#include "improc/improc.hpp"
#include <numeric>

ImProc::ImProc(ImCap *imCap, std::shared_ptr<logger> log)
    : conf(Config::conf), confPath(toml::get<std::string>(conf["improc"]["confPath"])),
      dataPath(toml::get<std::string>(conf["postproc"]["procDataPath"])), imCap(imCap),
      procData(new QueueFPS<std::vector<double>>(dataPath + "procDataQueue.txt")), lg(log) {
  lg->info("Initializing ImProc...");
  procData->out << "time (ms), ";
  procData->out << "yState1_1, yState1_2, yState1_3, yState2_1, yState2_2, yState2_3, ";
  procData->out << "y1, y2, y3, y4, y5, y6, ";
  procData->out << "yPrev1, yPrev2, yPrev3, yPrev4, yPrev5, yPrev6, ";

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
findInferredClstrs(const std::vector<std::vector<double>> &directClstrs,
                   std::shared_ptr<logger> lg) {
  int n = directClstrs.size();
  std::vector<std::vector<double>> newClstrs(n);

  for (int ch = 0; ch < n; ++ch) {
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
    lg->info("{} otherClstrs with valid yDirect selected for calculating yInferred in channel {}",
             m, ch);

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

void ImProc::initStates(int maxInitialDirectChs) {
  if (numInitDirectChs == maxInitialDirectChs)
    return;

  for (size_t ch = 0; ch < no; ++ch) {
    if (yDirect1[ch].has_value() && yDirect1[ch].value() < 0 &&
        numInitDirectChs < maxInitialDirectChs) {
      // Transition from inferred to direct
      yState1[ch] = true;
      yState2[ch] = true;
      numInitDirectChs++;
    }
  }
}

/*
 * update measurement given current state
 *   - direct = yDirect (if it's valid and close to previous measurement)
 *   - inferred = yInferred (if it's valid and close to previous measurement)
 */
void updateMeasCh(std::vector<double> &y, std::vector<double> &yPrev,
                  const std::vector<OptDouble> &yDirect, const std::vector<OptDouble> &yInferred,
                  std::vector<bool> &state) {
  int epsil = 56 / 2; // half of channel width

  for (size_t ch = 0; ch < y.size(); ++ch) {
    if (state[ch]) // Direct state
      y[ch] = (yDirect[ch].has_value() && std::abs(yDirect[ch].value() - yPrev[ch]) < epsil)
                  ? yDirect[ch].value()
                  : yPrev[ch];
    else // Inferred state
      y[ch] = (yInferred[ch].has_value() && std::abs(yInferred[ch].value() - yPrev[ch]) < epsil)
                  ? yInferred[ch].value()
                  : yPrev[ch];
  }
}

void ImProc::updateMeasAndStateOnZeroCross() {
  int epsil = 56 / 2; // half of channel width
  bool d2iTxRequested = false;

  for (int ch = 0; ch < no; ++ch) {
    if ((yState1[ch] && y1[ch] >= 0) || (yState2[ch] && y2[ch] >= 0)) {
      d2iTxRequested = true;
      lg->info("d2i transition requested on channel {}", ch);
    }
  }

  {
    // std::lock_guard<std::mutex> lock(zeroCrossMtx);
    zeroCross.assign(2 * no, false);
    bool y1Controlled = anyNonZeroR(0, no);
    bool y2Controlled = anyNonZeroR(no, 2 * no);
    if (d2iTxRequested && y1Controlled) {
      for (int ch = 0; ch < no; ++ch) {
        y2[ch] = yState2[ch] ? -1 : 1;
        zeroCross[no + ch] = true;
      }
    }
    if (d2iTxRequested && y2Controlled) {
      for (int ch = 0; ch < no; ++ch) {
        y1[ch] = yState1[ch] ? -1 : 1;
        zeroCross[ch] = true;
      }
    }
  }

  for (int ch = 0; ch < no; ++ch) {
    if (yState1[ch]) {
      if (y1[ch] >= 0)
        yState1[ch] = false;
    } else {
      OptDouble yDirect = argMinDist(directFgClstrs[ch], 0);
      if (yDirect.has_value() && yDirect.value() < 0 &&
          std::abs(yDirect.value() - y1[ch]) < epsil) {
        yState1[ch] = true;
        y1[ch] = yDirect.value();
      }
    }
    if (yState2[ch]) {
      if (y2[ch] >= 0)
        yState2[ch] = false;
    } else {
      OptDouble yDirect = argMinDist(directFgClstrs[ch], 0);
      if (yDirect.has_value() && yDirect.value() < 0 &&
          std::abs(yDirect.value() - y2[ch]) < epsil) {
        yState2[ch] = true;
        y2[ch] = yDirect.value();
      }
    }
  }
}

void ImProc::updateMeas() {
  // get closest direct & inferred clstr (yDirect and yInferred) to prev measurements if they exist
  yPrev1 = y1;
  yPrev2 = y2;
  for (int ch = 0; ch < no; ++ch) {
    directMeasAvail[ch] = !directFgClstrs[ch].empty();
    yDirect1[ch] = argMinDist(directFgClstrs[ch], yPrev1[ch]);
    yInferred1[ch] = argMinDist(inferredFgClstrs[ch], yPrev1[ch]);
    yDirect2[ch] = argMinDist(directFgClstrs[ch], yPrev2[ch]);
    yInferred2[ch] = argMinDist(inferredFgClstrs[ch], yPrev2[ch]);
    lg->info("ch: {}, yDirect1: {}, yInferred1: {}, yDirect2: {}, yInferred2: {}, yState1: {}, "
             "yState2 : {} ",
             ch, yDirect1[ch].value_or(1000), yInferred1[ch].value_or(1000),
             yDirect2[ch].value_or(1000), yInferred2[ch].value_or(1000), yState1[ch], yState2[ch]);
  }

  initStates(1);
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
    numInitDirectChs = 0;

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

      for (int ch = 0; ch < no; ++ch) {
        segAndOrientCh(tempFgMask, tempChFgMask, procFrameArr[ch], impConf.chROIs_[ch],
                       impConf.chWidth_);
        // find fgLocs in center column (excluding all +ve values except 2 rows above row Wch/2)
        cv::findNonZero(
            procFrameArr[ch](cv::Range(impConf.chWidth_ / 2 - 2, procFrameArr[ch].rows),
                             cv::Range(procFrameArr[ch].cols / 2, procFrameArr[ch].cols / 2 + 1)),
            fgLocs);
        // print fgLocs
        // for (int i = 0; i < fgLocs.size(); ++i)
        //   lg->info("ch: {}, i: {}, fgLocs: {}", ch, i, fgLocs[i]);

        findClusters(fgLocs, directFgClstrs[ch]);

        // coordinate transformation
        for (auto &clstrLoc : directFgClstrs[ch])
          clstrLoc = -clstrLoc + 2;
      }

      inferredFgClstrs = findInferredClstrs(directFgClstrs, lg);

      for (int ch = 0; ch < no; ++ch)
        for (int i = 0; i < directFgClstrs[ch].size(); ++i)
          lg->info("ch: {}, i: {}, directFgClstrs: {}", ch, i, directFgClstrs[ch][i]);
      for (int ch = 0; ch < no; ++ch)
        for (int i = 0; i < inferredFgClstrs[ch].size(); ++i)
          lg->info("ch: {}, i: {}, inferredFgClstrs: {}", ch, i, inferredFgClstrs[ch][i]);

      updateMeas();

      {
        std::lock_guard<std::mutex> lock(yMtx);
        y.clear();
        y.insert(y.end(), y1.begin(), y1.end());
        y.insert(y.end(), y2.begin(), y2.end());
      }
      // print y1 and y2
      for (int ch = 0; ch < no; ++ch)
        lg->info("ch: {}, y1: {}, y2: {}", ch, y1[ch], y2[ch]);

      procData->push(y);
      procFrameBuf.set(procFrameArr);
      procData->out << "\n";
    } catch (cv::Exception &e) {
      lg->error("Message: {}", e.what());
      lg->error("Type: {}", type_name<decltype(e)>());
    }
  }
}

bool ImProc::started() { return startedImProc; }

OptDouble ImProc::argMinDist(std::vector<double> &vec, double prevVal) {
  OptDouble argMin;
  if (vec.empty())
    return argMin; // empty opt indicates no argmin found (i.e. value is not up-to-date)

  return *std::min_element(vec.begin(), vec.end(), [prevVal](double a, double b) {
    return std::abs(a - prevVal) < std::abs(b - prevVal);
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
