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
    std::cout << "##For channel " << ch << ":\n";
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
    std::cout << m
              << " other channels (otherClstrs) with direct measurements selected for inference\n";

    // find all ways to choose an index from each of the m vectors in `otherClstrs`
    std::vector<int> indices(m, 0);
    findMeasCombs(otherClstrs, indices, 0, [&](const std::vector<int> &indices) {
      double sum = 0;
      for (int j = 0; j < m; ++j)
        sum += otherClstrs[j][indices[j]];
      newClstrs[ch].push_back(-sum / (n - 1));
    });

    // newClstrs[i].insert(newClstrs[i].end(), clstrs[i].begin(), clstrs[i].end());
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

void updateMeas(std::vector<double> &y, std::vector<double> &yPrev,
                const std::vector<bool> &directMeasAvail, const std::vector<double> &yDirect,
                const std::vector<double> &yInferred, std::vector<bool> &state, int &txCountDown,
                int &i2dOccurrences, bool &d2iTxOccurred) {
  yPrev = y;
  size_t d2iTxIdx = y.size(); // Initialize with an invalid index
  int maxInitI2D = 1;         // num of channels with direct measurements initially

  for (size_t ch = 0; ch < y.size(); ++ch) {
    if (state[ch]) { // Direct state
      y[ch] = std::min(yDirect[ch], 0.0);
      if (yDirect[ch] > 0 && txCountDown == 0) {
        // Transition from direct to inferred
        state[ch] = false;
        d2iTxOccurred = true;
        d2iTxIdx = ch;
        txCountDown = 20;
      }
    } else { // Inferred state
      y[ch] = yInferred[ch];
      if (directMeasAvail[ch] && yDirect[ch] < 0 && i2dOccurrences < maxInitI2D) {
        // Transition from inferred to direct
        state[ch] = true;
        i2dOccurrences++;
      }
    }
  }

  // If a direct->inferred transition occurred,
  // TODO all inferred states with yPrev and yDirect close to 0 become direct
  // (except the state that triggered the transition),
  // TODO all direct states close to 0 become inferred
  if (d2iTxOccurred) {
    for (size_t ch = 0; ch < y.size(); ++ch)
      if (!state[ch] && ch != d2iTxIdx)
        state[ch] = true;
    d2iTxOccurred = false;
    d2iTxIdx = y.size();
  }

  // Decrease transition countdown
  if (txCountDown > 0)
    --txCountDown;
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
    pBackSub =
        cv::createBackgroundSubtractorMOG2(impConf.bgSubHistory_, impConf.bgSubThres_, false);
    rectElement = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3));

    procFrameArr.assign(impConf.numChs_, cv::Mat());
    directFgClstrs.assign(impConf.numChs_, std::vector<double>());
    directMeasAvail.assign(impConf.numChs_, false);
    yState1.assign(impConf.numChs_, false);
    yState2.assign(impConf.numChs_, false);
    yDirect1.assign(impConf.numChs_, 0);
    yDirect2.assign(impConf.numChs_, 0);
    yInferred1.assign(impConf.numChs_, 0);
    yInferred2.assign(impConf.numChs_, 0);
    y1.assign(impConf.numChs_, 0);
    y2.assign(impConf.numChs_, 0);
    yPrev1.assign(impConf.numChs_, 0);
    yPrev2.assign(impConf.numChs_, 0);
    txCooldown1 = 0, txCooldown2 = 0, i2dOccurrences1 = 0, i2dOccurrences2 = 0;
    d2iTxOccurred1 = false, d2iTxOccurred2 = false;

    // TODO handle yMax initialization for non-90-degree channels
    yMax.clear();
    for (int ch = 0; ch < impConf.numChs_; ch++)
      (ch == 0) ? yMax.push_back(impConf.getChROIs()[ch].height - impConf.getChWidth())
                : yMax.push_back(impConf.getChROIs()[ch].width - impConf.getChWidth());
    procThread = std::thread(&ImProc::start, this);
    procThread.detach();
  }
}

void ImProc::stopThread() {
  if (started()) {
    lg->info("Stopping image processing...");
    startedImProc = false;
  }
}

void ImProc::start() {
  while (started()) {
    preFrame = imCap->getFrame();
    try {
      // Timer t("ImProc");

      // update foreground mask based on background threshold,
      // (optionally update background model at preset learning rate),
      // then apply dilation/erosion to remove noise
      pBackSub->apply(preFrame, fgMask, 0);
      cv::dilate(fgMask, tempFgMask, rectElement, cv::Point(-1, -1), 2);
      cv::erode(tempFgMask, tempFgMask, rectElement, cv::Point(-1, -1), 4);
      cv::dilate(tempFgMask, tempFgMask, rectElement, cv::Point(-1, -1), 2);
      skeletonize(tempFgMask, rectElement);
      cv::dilate(tempFgMask, tempFgMask, rectElement);

      for (int ch = 0; ch < impConf.numChs_; ++ch) {
        segAndOrientCh(tempFgMask, tempChFgMask, procFrameArr[ch], impConf.chROIs_[ch],
                       impConf.chWidth_);
        // find fgLocs in center column
        cv::findNonZero(procFrameArr[ch].col(procFrameArr[ch].cols / 2), fgLocs);
        // print fgLocs
        // for (int i = 0; i < fgLocs.size(); ++i)
        //   info("ch: {}, i: {}, fgLocs: {}", ch, i, fgLocs[i]);

        findClusters(fgLocs, directFgClstrs[ch]);

        // coordinate transformation
        for (auto &clstrLoc : directFgClstrs[ch])
          clstrLoc = -(clstrLoc - impConf.chWidth_ / 2.0);
      }

      inferredFgClstrs = findInferredClstrs(directFgClstrs);
      // print each cluster
      // info("directFgClstrs size: {}", directFgClstrs.size());
      for (int ch = 0; ch < impConf.numChs_; ++ch) {
        // print size of each cluster
        // info("ch: {}, directFgClstrs size: {}", ch, directFgClstrs[ch].size());
        for (int i = 0; i < directFgClstrs[ch].size(); ++i)
          lg->info("ch: {}, i: {}, directFgClstrs: {}", ch, i, directFgClstrs[ch][i]);
      }

      // for (int ch = 0; ch < impConf.numChs_; ++ch)
      //   for (int i = 0; i < directFgClstrs[ch].size(); ++i)
      //     info("ch: {}, i: {}, inferredFgClstrs: {}", ch, i, inferredFgClstrs[ch][i]);

      for (int ch = 0; ch < impConf.numChs_; ++ch) {
        directMeasAvail[ch] = !directFgClstrs[ch].empty();
        yDirect1[ch] = minDist(directFgClstrs[ch], yPrev1[ch]);
        yInferred1[ch] = minDist(inferredFgClstrs[ch], yPrev1[ch]);
        yDirect2[ch] = minDist(directFgClstrs[ch], yPrev2[ch]);
        yInferred2[ch] = minDist(inferredFgClstrs[ch], yPrev2[ch]);
      }
      updateMeas(y1, yPrev1, directMeasAvail, yDirect1, yInferred1, yState1, txCooldown1,
                 i2dOccurrences1, d2iTxOccurred1);
      updateMeas(y2, yPrev2, directMeasAvail, yDirect2, yInferred2, yState2, txCooldown2,
                 i2dOccurrences2, d2iTxOccurred2);

      // for (int ch = 0; ch < impConf.numChs_; ++ch) {
      //   double y1Direct = minDist(directFgClstrs[ch], yPrev1[ch]);
      //   double y1Inferred = minDist(inferredFgClstrs[ch], yPrev1[ch]);
      //   double y2Direct = minDist(directFgClstrs[ch], yPrev2[ch]);
      //   double y2Inferred = minDist(inferredFgClstrs[ch], yPrev2[ch]);

      //   lg->info("ch: {}, y1Direct: {}, y1Inferred: {}, y2Direct: {}, y2Inferred: {}", ch,
      //   y1Direct,
      //            y1Inferred, y2Direct, y2Inferred);

      //   // zero-crossing detection
      //   if ((!directFgClstrs[ch].empty()) && (yPrev1[ch] < 0)) { // direct
      //     yState1[ch] = true;
      //     y1[ch] = y1Direct;
      //   } else { // inferred
      //     yState1[ch] = false;
      //     y1[ch] = ((!directFgClstrs[ch].empty()) && (y1Direct < 0)) ? y1Direct : y1Inferred;
      //   }
      //   if ((!directFgClstrs[ch].empty()) && (yPrev2[ch] < 0)) { // direct
      //     yState2[ch] = true;
      //     y2[ch] = y2Direct;
      //   } else { // inferred
      //     yState2[ch] = false;
      //     y2[ch] = ((!directFgClstrs[ch].empty()) && (y2Direct < 0)) ? y2Direct : y2Inferred;
      //   }

      //   yPrev1[ch] = y1[ch];
      //   yPrev2[ch] = y2[ch];
      // }

      /*
      ** - each channel can be in 2 states: direct or inferred
      ** - these states correspond to the 2 types of measurements
      ** - direct measurements are always -ve, inferred are always +ve
      ** - when a direct measurement crosses from -ve to +ve,
      **     it becomes an inferred measurement and we discard any future direct measurements
      ** - when an inferred measurement crosses from +ve to -ve,
      **     it becomes a direct measurement and we discard any future inferred measurements
      */
      rstOnZeroCross();

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
    if (r[i] == 0)
      return false;
  return true;
}

bool ImProc::anyZeroCross(const std::vector<double> &vec1, const std::vector<double> &vec2) {
  if (vec1.size() != vec2.size())
    throw std::invalid_argument("Vectors are not the same size");

  auto it = std::mismatch(vec1.begin(), vec1.end(), vec2.begin(),
                          [this](double val1, double val2) { return (val1 >= 0) == (val2 >= 0); });

  return it.first != vec1.end();
}

void ImProc::rstOnZeroCross() {
  // if y1 is being controlled and any y2 crosses zero, reset y2
  if (anyNonZeroR(0, impConf.numChs_))
    if (anyZeroCross(y2, yPrev2))
      std::fill(y2.begin(), y2.end(), 0);
  // if y2 is being controlled and any y1 crosses zero, reset y1
  if (anyNonZeroR(impConf.numChs_, 2 * impConf.numChs_))
    if (anyZeroCross(y1, yPrev1))
      std::fill(y1.begin(), y1.end(), 0);
}

void ImProc::setR(double currTraj[2 * MAX_NO]) {
  for (int i = 0; i < 2 * impConf.numChs_; ++i)
    r[i] = currTraj[i];
}

std::vector<double> ImProc::getY() {
  std::lock_guard<std::mutex> lock(yMtx);
  return y;
}

cv::Mat ImProc::getProcFrame(int idx) { return procFrameBuf.get()[idx]; }

void ImProc::clearData() { procData->clear(); }
