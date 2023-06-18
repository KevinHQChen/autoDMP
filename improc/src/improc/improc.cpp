#include "improc/improc.hpp"
#include <numeric>

void combinations(std::vector<std::vector<double>> &clstrs, std::vector<std::vector<double>> &yVecs,
                  std::vector<double> y, int ch, int selChs) {
  // selChs base case: we've selected enough chs,
  // fill remaining chs with -1 and add y to yVecs
  if (selChs == 0) {
    while (y.size() < clstrs.size()) {
      y.push_back(-1);
    }
    yVecs.push_back(y);
    return;
  }

  // ch base case: we've processed all chs, don't add this combination to yVecs
  if (ch == clstrs.size())
    return;

  // recursive case: call combinations for each cluster in current ch
  for (int i = 0; i < clstrs[ch].size(); ++i) {
    // Add current cluster location to y
    y.push_back(clstrs[ch][i]);
    // Recursively call combinations with the next ch and one less selChs
    combinations(clstrs, yVecs, y, ch + 1, selChs - 1);
    // reset y for next iteration
    y.pop_back();
  }

  // recursive case: call combinations for next ch, assume no clusters exist in current ch
  y.push_back(-1);
  combinations(clstrs, yVecs, y, ch + 1, selChs);
  y.pop_back();
}

std::vector<std::vector<double>> findCombinations(std::vector<std::vector<double>> &clstrs) {
  std::vector<std::vector<double>> yVecs;
  // Call combinations for each possible selChs up to numChs - 1
  for (int selChs = 1; selChs < clstrs.size(); ++selChs) {
    std::vector<double> y;
    combinations(clstrs, yVecs, y, 0, selChs);
  }
  // return yVecs which now contains all combinations
  return yVecs;
}

void applyKCL(std::vector<std::vector<double>> &yVecs) {
  // Apply KCl to approximate unmeasured channels in yVecs
  for (auto &y : yVecs) {
    int sum = 0;

    for (auto &num : y)
      if (num >= 0)
        sum += num;

    for (auto &num : y) {
      if (num >= 0)
        num = -num;
      else
        num = (y.size() - 1) == 0 ? 0 : sum / (y.size() - 1.0);
    }
  }
}

double l2Norm(const std::vector<double> &a, const std::vector<double> &b) {
  return std::sqrt(std::inner_product(a.begin(), a.end(), b.begin(), 0.0, std::plus<double>(),
                                      [](double x, double y) { return (x - y) * (x - y); }));
}

std::vector<double> minL2Norm(const std::vector<std::vector<double>> &yVecs,
                              const std::vector<double> &yPrev) {
  auto yMinL2 =
      std::min_element(yVecs.begin(), yVecs.end(),
                       [&yPrev](const std::vector<double> &y1, const std::vector<double> &y2) {
                         return l2Norm(y1, yPrev) < l2Norm(y2, yPrev);
                       });
  return *yMinL2;
}

ImProc::ImProc(ImCap *imCap)
    : conf(Config::conf), confPath(toml::get<std::string>(conf["improc"]["confPath"])),
      dataPath(toml::get<std::string>(conf["postproc"]["procDataPath"])), imCap(imCap),
      procData(new QueueFPS<std::vector<double>>(dataPath + "procDataQueue.txt")) {
  for (int ch = 0; ch < impConf.numChs_; ch++)
    procData->out << "time (ms), maxLoc.x (px), maxLoc.y (px)\n";

  // save images with proper format PNG, CV_16UC1
  compParams.push_back(cv::IMWRITE_PNG_COMPRESSION);
  compParams.push_back(0);
}

ImProc::~ImProc() {
  stopProcThread();
  delete procData;
}

// load channel/template images, bounding boxes from file
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

void ImProc::startProcThread() {
  if (!started()) {
    info("Starting image processing...");
    startedImProc = true;
    pBackSub =
        cv::createBackgroundSubtractorMOG2(impConf.bgSubHistory_, impConf.bgSubThres_, false);
    for (int ch = 0; ch < impConf.numChs_; ch++) {
      procFrameArr.push_back(cv::Mat());
      fgClstrs.push_back(std::vector<double>());
      y1.push_back(0);
      y2.push_back(0);
      yPrev1.push_back(0);
      yPrev2.push_back(0);
      // TODO add support for non-90-degree channels
      if (ch == 0)
        yMax.push_back(impConf.getChROIs()[ch].height - impConf.getChWidth());
      else
        yMax.push_back(impConf.getChROIs()[ch].width - impConf.getChWidth());
    }
    procThread = std::thread(&ImProc::start, this);
    procThread.detach();
  }
}

void ImProc::stopProcThread() {
  if (started()) {
    info("Stopping image processing...");
    startedImProc = false;
    procFrameArr.clear();
    fgClstrs.clear();
    y1.clear();
    y2.clear();
    yPrev1.clear();
    yPrev2.clear();
    yMax.clear();
    if (procThread.joinable())
      procThread.join();
  }
}

void ImProc::start() {
  while (started()) {
    preFrame = imCap->getFrame();
    try {
      auto startTime = high_resolution_clock::now();
      auto stopTime = high_resolution_clock::now();

      // update foreground mask based on background threshold,
      // (optionally update background model at preset learning rate),
      // then apply dilation/erosion to remove noise
      pBackSub->apply(preFrame, fgMask, 0);
      cv::dilate(fgMask, tempFgMask, cv::Mat(), cv::Point(-1, -1), 2);
      cv::erode(tempFgMask, tempFgMask, cv::Mat(), cv::Point(-1, -1), 6);
      cv::dilate(tempFgMask, tempFgMask, cv::Mat(), cv::Point(-1, -1), 2);

      for (int ch = 0; ch < impConf.numChs_; ++ch) {
        segAndOrientCh(tempFgMask, procFrameArr[ch], impConf.chROIs_[ch], impConf.chWidth_);
        cv::findNonZero(
            procFrameArr[ch](cv::Rect(0, yMax[ch], impConf.getChWidth(), impConf.getChWidth())),
            fgLocs);
        findClusters(fgLocs, fgClstrs[ch], 0);
      }

      yVecs = findCombinations(fgClstrs);
      applyKCL(yVecs);
      y1 = minL2Norm(yVecs, yPrev1);
      y2 = minL2Norm(yVecs, yPrev2);
      yPrev1 = y1;
      yPrev2 = y2;

      y.clear();
      y.insert( y.end(), y1.begin(), y1.end() );
      y.insert( y.end(), y2.begin(), y2.end() );
      procData->push(y);
      procFrameBuf.set(procFrameArr);

      stopTime = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stopTime - startTime);
      // info("imProc duration: {}", duration.count());
    } catch (cv::Exception &e) {
      error("Message: {}", e.what());
      error("Type: {}", type_name<decltype(e)>());
    }
  }
}

void ImProc::segAndOrientCh(cv::Mat &srcImg, cv::Mat &destImg, RotRect &chROI, int &chWidth) {
  cv::Mat tmpImg = srcImg(chROI);
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

void ImProc::findClusters(const std::vector<cv::Point> &fgLocs, std::vector<double> &clusters,
                          int tolerance) {
  // we assume fgLocs is sorted by y-coordinate, regardless whether it's ascending or descending
  // (just has to be consistent), but this is undocumented

  clusters.clear();
  std::optional<int> start;
  std::optional<int> lastFgIdx;

  if (fgLocs.empty()) {
    clusters = {};
    return;
  }

  for (const auto &pt : fgLocs) {
    if (!start)
      start = pt.y;
    else if (lastFgIdx && pt.y - *lastFgIdx > tolerance) {
      clusters.push_back((*start + *lastFgIdx) / 2.0);
      start = pt.y;
    }
    lastFgIdx = pt.y;
  }

  // Handle the case where the last cluster is at the end of the image
  if (start && lastFgIdx && *lastFgIdx - *start <= tolerance)
    clusters.push_back((*start + *lastFgIdx) / 2.0);
}

bool ImProc::started() { return startedImProc; }

cv::Mat ImProc::getProcFrame(int idx) { return procFrameBuf.get()[idx]; }

void ImProc::clearProcData() { procData->clear(); }
