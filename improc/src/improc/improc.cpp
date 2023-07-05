#include "improc/improc.hpp"
#include <numeric>

ImProc::ImProc(ImCap *imCap, std::shared_ptr<logger> log)
    : conf(Config::conf), confPath(toml::get<std::string>(conf["improc"]["confPath"])),
      dataPath(toml::get<std::string>(conf["postproc"]["procDataPath"])), imCap(imCap),
      procData(new QueueFPS<std::vector<double>>(dataPath + "procDataQueue.txt")), lg(log) {
  lg->info("Initializing ImProc...");
  // for (int ch = 0; ch < impConf.numChs_; ++ch)
  //   procData->out << "time (ms), fgLoc.y (px)\n";

  // save images with proper format PNG, CV_16UC1
  compParams.push_back(cv::IMWRITE_PNG_COMPRESSION);
  compParams.push_back(0);
}

ImProc::~ImProc() {
  lg->info("Terminating ImProc...");
  stopThread();
  delete procData;
}

void findMeasCombs(const std::vector<std::vector<double>> &sets, std::vector<int> &indices,
                   int setIndex, std::function<void(const std::vector<int> &)> callback) {
  if (setIndex == sets.size()) {
    callback(indices);
    return;
  }
  for (int i = 0; i < sets[setIndex].size(); ++i) {
    indices[setIndex] = i;
    findMeasCombs(sets, indices, setIndex + 1, callback);
  }
}

void findChCombs(int offset, int j, int m, std::vector<int> &combination,
                 std::vector<std::vector<int>> &combinations) {
  if (j == 0) {
    combinations.push_back(combination);
    return;
  }
  for (int i = offset; i <= m - j; ++i) {
    combination.push_back(i);
    findChCombs(i + 1, j - 1, m, combination, combinations);
    combination.pop_back();
  }
}

std::vector<std::vector<double>>
findInferredClstrs(const std::vector<std::vector<double>> &directClstrs) {
  int n = directClstrs.size();
  std::vector<std::vector<double>> newClstrs(n);

  for (int ch = 0; ch < n; ++ch) {
    // std::cout << "##For channel " << ch << ":\n";
    // make a copy of directClstrs_,
    // remove i-th channel (we only use other channels to infer clstr locations in current channel)
    // remove empty channels (i.e. channels w/o direct clstrs),
    // remove +ve clstrs in each channel (we won't use them for inference),
    std::vector<std::vector<double>> otherClstrs(directClstrs);
    otherClstrs.erase(otherClstrs.begin() + ch);
    for (auto &ch : otherClstrs)
      ch.erase(std::remove_if(ch.begin(), ch.end(), [](double x) { return x > 0; }), ch.end());
    otherClstrs.erase(std::remove_if(otherClstrs.begin(), otherClstrs.end(),
                                     [](const std::vector<double> &v) { return v.empty(); }),
                      otherClstrs.end());

    int m = otherClstrs.size();
    // std::cout << m
    //           << " other channels (otherClstrs) with direct measurements selected for
    //           inference\n";

    // find all ways to choose an index from each of the m vectors in `otherClstrs`
    std::vector<int> indices(m, 0);
    findMeasCombs(otherClstrs, indices, 0, [&](const std::vector<int> &indices) {
      double sum = 0;
      for (int j = 0; j < m; ++j)
        sum += otherClstrs[j][indices[j]];
      newClstrs[ch].push_back(-sum / (n - 1));
    });
  }

  return newClstrs;
}

void skeletonize(cv::Mat &img, cv::Mat element) {
  // skeletonize
  if (cv::countNonZero(img) == img.total())
    img.setTo(cv::Scalar(0));
  cv::Mat skel(img.size(), CV_8UC1, cv::Scalar(0));
  cv::Mat temp;
  cv::Mat eroded;
  bool done;
  do {
    cv::erode(img, eroded, element);
    cv::dilate(eroded, temp, element); // temp = open(img)
    cv::subtract(img, temp, temp);
    cv::bitwise_or(skel, temp, skel);
    eroded.copyTo(img);
    done = (cv::countNonZero(img) == 0);
  } while (!done);
  img = skel;
}

bool ImProc::initStates() {
  if (doneInit)
    return true;

  int maxInitialDirectChs = no - 2;

  for (size_t ch = 0; ch < no; ++ch) {
    if (directMeasAvail[ch] && (yDirect1[ch] < 0) && (numInitDirectChs < maxInitialDirectChs)) {
      // Transition from inferred to direct
      yState1[ch] = true;
      yState2[ch] = true;
      numInitDirectChs++;
    }
  }

  if (numInitDirectChs == maxInitialDirectChs) {
    doneInit = true;
    return true;
  }

  return false;
}

/*
 * update measurement given current state
 *   - direct = yDirect (or yPrev if yDirect is too far away)
 *   - inferred = yInferred
 */
void updateMeasCh(std::vector<double> &y, std::vector<double> &yPrev,
                  const std::vector<double> &yDirect, const std::vector<double> &yInferred,
                  std::vector<bool> &state) {
  int epsilon = 56 / 2; // half of channel width

  for (size_t ch = 0; ch < y.size(); ++ch) {
    if (state[ch]) // Direct state
      y[ch] = ((yDirect[ch] > 0) && (yPrev[ch] < 0) && (yDirect[ch] - yPrev[ch] > epsilon))
                  ? yPrev[ch]
                  : yDirect[ch];
    else // Inferred state
      y[ch] = yInferred[ch];
  }
}

void ImProc::updateMeasAndStateOnZeroCross() {
  int epsil = 56 / 2; // half of channel width
  bool d2iTxRequested = false;
  size_t d2iCh = no;
  std::vector<bool> yStatePrev1(yState1);
  std::vector<bool> yStatePrev2(yState2);

  for (int ch = 0; ch < no; ++ch) {
    if (((yState1[ch] && y1[ch] >= 0) || (yState2[ch] && y2[ch] >= 0)) && (txCooldown == 0)) {
      d2iTxRequested = true;
      d2iCh = ch;
      txCooldown = 20;
    }
  }
  if (txCooldown > 0)
    --txCooldown;

  // {
  // std::lock_guard<std::mutex> lock(zeroCrossMtx);
  zeroCross.assign(2 * no, false);
  bool y1Controlled = anyNonZeroR(0, no);
  bool y2Controlled = anyNonZeroR(no, 2 * no);
  if (d2iTxRequested && y1Controlled) {
    for (int ch = 0; ch < no; ++ch) {
      yState2[ch] = (ch == d2iCh) ? false : true;
      yDirect2[ch] = minDist(directFgClstrs[ch], 0);
      yInferred2[ch] = minDist(inferredFgClstrs[ch], 0);
      yPrev2[ch] = yInferred2[ch];
      zeroCross[no + ch] = true;
    }
    updateMeasCh(y2, yPrev2, yDirect2, yInferred2, yState2);
  }
  if (d2iTxRequested && y2Controlled) {
    for (int ch = 0; ch < no; ++ch) {
      yState1[ch] = (ch == d2iCh) ? false : true;
      yDirect1[ch] = minDist(directFgClstrs[ch], 0);
      yInferred1[ch] = minDist(inferredFgClstrs[ch], 0);
      yPrev1[ch] = yInferred1[ch];
      zeroCross[ch] = true;
    }
    updateMeasCh(y1, yPrev1, yDirect1, yInferred1, yState1);
  }
  // }

  if (d2iTxRequested && !y2Controlled) {
    for (size_t ch = 0; ch < no; ++ch) {
      if (yState1[ch] && (y1[ch] >= 0 || std::abs(y1[ch]) < epsil))
        yState1[ch] = false;
      else if (!yState1[ch] && (std::abs(y1[ch]) < epsil) && (std::abs(yDirect1[ch]) < epsil))
        yState1[ch] = true;
    }
  }
  if (d2iTxRequested && !y1Controlled) {
    for (size_t ch = 0; ch < no; ++ch) {
      if (yState2[ch] && (y2[ch] >= 0 || std::abs(y2[ch]) < epsil))
        yState2[ch] = false;
      else if (!yState2[ch] && (std::abs(y2[ch]) < epsil) && (std::abs(yDirect2[ch]) < epsil))
        yState2[ch] = true;
    }
  }
}

void ImProc::updateMeas() {
  /*
   * get closest direct & inferred clstr (yDirect and yInferred)
   * to previous measurements,
   * if any of them don't exist, just use the previous measurement
   */
  yPrev1 = y1;
  yPrev2 = y2;
  for (int ch = 0; ch < no; ++ch) {
    directMeasAvail[ch] = !directFgClstrs[ch].empty();
    yDirect1[ch] = minDist(directFgClstrs[ch], yPrev1[ch]);
    yInferred1[ch] = minDist(inferredFgClstrs[ch], yPrev1[ch]);
    yDirect2[ch] = minDist(directFgClstrs[ch], yPrev2[ch]);
    yInferred2[ch] = minDist(inferredFgClstrs[ch], yPrev2[ch]);
    // lg->info("ch: {}, yDirect1: {}, yInferred1: {}, yDirect2: {}, yInferred2: {}, yState1:
    // {},
    // "
    //          "yState2: {}",
    //          ch, yDirect1[ch], yInferred1[ch], yDirect2[ch], yInferred2[ch], yState1[ch],
    //          yState2[ch]);
  }

  initStates();
  updateMeasCh(y1, yPrev1, yDirect1, yInferred1, yState1);
  updateMeasCh(y2, yPrev2, yDirect2, yInferred2, yState2);
  updateMeasAndStateOnZeroCross();
}

void ImProc::loadConfig() {
  // load bboxes from file
  ordered_value v = toml::parse<toml::discard_comments, tsl::ordered_map>(confPath + "config.toml");
  impConf.from_toml(v);
}

void ImProc::saveConfig() {
  // save bboxes to file
  ordered_value newImProcConfig(impConf);
  std::ofstream out(confPath + "config.toml");
  out << toml::format(newImProcConfig);
  out.close();
}

void ImProc::startThread() {
  if (!started()) {
    lg->info("Starting image processing...");
    startedImProc = true;
    no = impConf.numChs_;
    pBackSub =
        cv::createBackgroundSubtractorMOG2(impConf.bgSubHistory_, impConf.bgSubThres_, false);
    rectElement = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3));

    procFrameArr.assign(no, cv::Mat());
    directFgClstrs.assign(no, std::vector<double>());
    y1.assign(no, 0);
    y2.assign(no, 0);
    yPrev1.assign(no, 0);
    yPrev2.assign(no, 0);
    yState1.assign(no, false);
    yState2.assign(no, false);
    yDirect1.assign(no, 0);
    yDirect2.assign(no, 0);
    yInferred1.assign(no, 0);
    yInferred2.assign(no, 0);
    directMeasAvail.assign(no, false);
    // TODO handle yMax initialization for non-90-degree channels
    yMax.clear();
    for (int ch = 0; ch < no; ch++)
      (ch == 0) ? yMax.push_back(impConf.getChROIs()[ch].height - impConf.getChWidth())
                : yMax.push_back(impConf.getChROIs()[ch].width - impConf.getChWidth());

    zeroCross.assign(2 * no, false);
    txCooldown = 0, numInitDirectChs = 0;
    doneInit = false;

    procThread = std::thread(&ImProc::start, this);
    procThread.detach();
  }
}

void ImProc::stopThread() {
  if (started()) {
    lg->info("Stopping image processing...");
    startedImProc = false;
    std::lock_guard<std::mutex> guard(imProcMtx); // wait for thread to finish
  }
}

void ImProc::start() {
  while (started()) {
    std::lock_guard<std::mutex> guard(imProcMtx);
    preFrame = imCap->getFrame();
    try {
      // Timer t("ImProc");

      // update foreground mask based on background threshold,
      // (optionally update background model at preset learning rate),
      // then apply dilation/erosion to remove noise
      pBackSub->apply(preFrame, fgMask, 0);
      cv::dilate(fgMask, tempFgMask, rectElement, cv::Point(-1, -1), 2);
      cv::erode(tempFgMask, tempFgMask, rectElement, cv::Point(-1, -1), 6);
      cv::dilate(tempFgMask, tempFgMask, rectElement, cv::Point(-1, -1), 2);
      skeletonize(tempFgMask, rectElement);
      cv::dilate(tempFgMask, tempFgMask, rectElement);

      for (int ch = 0; ch < impConf.numChs_; ++ch) {
        segAndOrientCh(tempFgMask, tempChFgMask, procFrameArr[ch], impConf.chROIs_[ch],
                       impConf.chWidth_);
        // find fgLocs in center column (excluding all +ve values except first 2 rows)
        cv::findNonZero(
            procFrameArr[ch](cv::Range(impConf.chWidth_ / 2 - 2, procFrameArr[ch].rows),
                             cv::Range(procFrameArr[ch].cols / 2, procFrameArr[ch].cols / 2 + 1)),
            fgLocs);
        // print fgLocs
        // for (int i = 0; i < fgLocs.size(); ++i)
        //   info("ch: {}, i: {}, fgLocs: {}", ch, i, fgLocs[i]);

        findClusters(fgLocs, directFgClstrs[ch]);

        // coordinate transformation
        for (auto &clstrLoc : directFgClstrs[ch])
          clstrLoc = -clstrLoc + 2;
      }

      inferredFgClstrs = findInferredClstrs(directFgClstrs);
      // print each cluster
      // info("directFgClstrs size: {}", directFgClstrs.size());
      // for (int ch = 0; ch < impConf.numChs_; ++ch) {
      //   // print size of each cluster
      //   // info("ch: {}, directFgClstrs size: {}", ch, directFgClstrs[ch].size());
      //   for (int i = 0; i < directFgClstrs[ch].size(); ++i)
      //     lg->info("ch: {}, i: {}, directFgClstrs: {}", ch, i, directFgClstrs[ch][i]);
      // }

      // for (int ch = 0; ch < impConf.numChs_; ++ch)
      //   for (int i = 0; i < directFgClstrs[ch].size(); ++i)
      //     info("ch: {}, i: {}, inferredFgClstrs: {}", ch, i, inferredFgClstrs[ch][i]);

      updateMeas();

      // if (zeroCross1)
      //   lg->info("Zero crossing occurred in y1");
      // if (zeroCross2)
      //   lg->info("Zero crossing occurred in y2");
      // for (int ch = 0; ch < impConf.numChs_; ++ch) {
      //   lg->info("yState1: {}", yState1[ch]);
      //   lg->info("yState2: {}", yState2[ch]);
      // }

      {
        std::lock_guard<std::mutex> lock(yMtx);
        y.clear();
        y.insert(y.end(), y1.begin(), y1.end());
        y.insert(y.end(), y2.begin(), y2.end());
      }

      procData->push(y);
      procFrameBuf.set(procFrameArr);
      // print y1 and y2
      // for (int ch = 0; ch < impConf.numChs_; ++ch)
      //   lg->info("ch: {}, y1: {}, y2: {}", ch, y1[ch], y2[ch]);
    } catch (cv::Exception &e) {
      lg->error("Message: {}", e.what());
      lg->error("Type: {}", type_name<decltype(e)>());
    }
  }
}

bool ImProc::started() { return startedImProc; }

double ImProc::minDist(std::vector<double> &vec, double value) {
  if (vec.empty())
    return value;

  return *std::min_element(vec.begin(), vec.end(), [value](double a, double b) {
    return std::abs(a - value) < std::abs(b - value);
  });
}

void ImProc::segAndOrientCh(cv::Mat &srcImg, cv::Mat &tmpImg, cv::Mat &destImg, RotRect &chROI,
                            int &chWidth) {
  tmpImg = srcImg(chROI);
  switch (chROI.angle) {
  case 0:
    destImg = tmpImg;
    break;
  case 90:
    cv::rotate(tmpImg, destImg, cv::ROTATE_90_COUNTERCLOCKWISE);
    break;
  case -90:
    cv::rotate(tmpImg, destImg, cv::ROTATE_90_CLOCKWISE);
    break;
  case 180:
    cv::rotate(tmpImg, destImg, cv::ROTATE_180);
    break;
  default:
    rotateMat(tmpImg, destImg, chROI.angle);
    int x = destImg.cols / 2 - chWidth / 2;
    int y = 0;
    int chHeight = destImg.rows;
    destImg = destImg(cv::Rect(x, y, chWidth, chHeight));
    break;
  }
}

void ImProc::findClusters(const std::vector<cv::Point> &fgLocs, std::vector<double> &clusters) {
  clusters.clear();
  if (fgLocs.empty())
    return;

  int clusterStart = fgLocs[0].y;
  int clusterEnd = fgLocs[0].y;
  for (size_t i = 1; i < fgLocs.size(); ++i) {
    if (fgLocs[i].y == clusterEnd + 1) // This point is part of the current cluster.
      clusterEnd = fgLocs[i].y;
    else {
      // This point starts a new cluster. Compute the center of the current cluster.
      double clusterCenter = (clusterStart + clusterEnd) / 2.0;
      clusters.push_back(clusterCenter);

      // Start a new cluster.
      clusterStart = fgLocs[i].y;
      clusterEnd = fgLocs[i].y;
    }
  }

  // Compute the center of the last cluster.
  double clusterCenter = (clusterStart + clusterEnd) / 2.0;
  clusters.push_back(clusterCenter);
}

bool ImProc::anyNonZeroR(std::size_t start, std::size_t end) {
  for (std::size_t i = start; i < end; ++i)
    if (r[i] != 0)
      return true;
  return false;
}

void ImProc::setR(double currTraj[2 * MAX_NO]) {
  for (int i = 0; i < 2 * no; ++i)
    r[i] = currTraj[i];
}

std::vector<double> ImProc::getY() {
  std::lock_guard<std::mutex> lock(yMtx);
  return y;
}

unsigned char ImProc::getZeroCross(int ch) {
  std::lock_guard<std::mutex> lock(zeroCrossMtx);
  return zeroCross[ch];
}

cv::Mat ImProc::getProcFrame(int idx) { return procFrameBuf.get()[idx]; }

void ImProc::clearData() { procData->clear(); }
