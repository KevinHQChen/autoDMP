#include "improc/improc.hpp"
#include <numeric>

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

double minDist(std::vector<double> &vec, double value) {
  if (vec.empty())
    return value;

  return *std::min_element(vec.begin(), vec.end(), [value](double a, double b) {
    return std::abs(a - value) < std::abs(b - value);
  });
}

std::vector<std::vector<double>> findOtherClstrs(const std::vector<std::vector<double>> &clstrs) {
  int n = clstrs.size();
  std::vector<std::vector<double>> newClstrs(n);

  for (int i = 0; i < n; ++i) {
    // std::cout << "##For clstrs " << i << ":\n";
    // store all non-empty elements of `newClstrs` except `newClstrs[i]` in `otherClstrs`
    std::vector<std::vector<double>> otherClstrs(clstrs);
    otherClstrs.erase(otherClstrs.begin() + i);
    otherClstrs.erase(std::remove_if(otherClstrs.begin(), otherClstrs.end(),
                                     [](const std::vector<double> &v) { return v.empty(); }),
                      otherClstrs.end());
    int m = otherClstrs.size();
    // std::cout << m << " otherClstrs selected\n";

    for (int j = 1; j < m + 1; ++j) {
      // std::cout << j << " selected channels\n";
      // find all distinct combinations of j indices from all indices of `otherClstrs`
      std::vector<std::vector<int>> clstrIndices;
      std::vector<int> clstrIdx;
      findChCombs(0, j, m, clstrIdx, clstrIndices);
      int l = clstrIndices.size();
      // std::cout << l << " set combinations found\n";
      // for (int i = 0; i < l; ++i) {
      //   std::cout << "clstrIndices[" << i << "] = {";
      //   for (int j = 0; j < clstrIndices[i].size(); ++j) {
      //     std::cout << clstrIndices[i][j];
      //     if (j != clstrIndices[i].size() - 1)
      //       std::cout << ", ";
      //   }
      //   std::cout << "}\n";
      // }

      for (int k = 0; k < l; ++k) {
        // std::cout << "combination " << k << ", ";
        // given the k-th combination of j indices, store the elements of `otherClstrs` at the
        // corresponding indices in `selClstrs`
        std::vector<std::vector<double>> selClstrs;
        for (int x = 0; x < j; ++x) {
          // std::cout << "selClstrs idx " << clstrIndices[k][x];
          selClstrs.push_back(otherClstrs[clstrIndices[k][x]]);
          // for (int i = 0; i < otherClstrs[clstrIndices[k][x]].size(); ++i)
          //   std::cout << ", val " << otherClstrs[clstrIndices[k][x]][i];
          // std::cout << "\n";
        }

        // find all distinct combinations of j elements from `selClstrs`
        std::vector<int> indices(j, 0);
        findMeasCombs(selClstrs, indices, 0, [&](const std::vector<int> &indices) {
          double sum = 0;
          for (int j_ = 0; j_ < j; ++j_)
            sum += selClstrs[j_][indices[j_]];
          newClstrs[i].push_back(-sum / (n - 1));
        });
      }
    }

    newClstrs[i].insert(newClstrs[i].end(), clstrs[i].begin(), clstrs[i].end());
  }

  return newClstrs;
}

ImProc::ImProc(ImCap *imCap)
    : conf(Config::conf), confPath(toml::get<std::string>(conf["improc"]["confPath"])),
      dataPath(toml::get<std::string>(conf["postproc"]["procDataPath"])), imCap(imCap),
      currEv(nullEv), procData(new QueueFPS<std::vector<double>>(dataPath + "procDataQueue.txt")) {
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
    std::lock_guard<std::mutex> guard(dataMtx); // wait for procThread to finish
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
    std::lock_guard<std::mutex> guard(dataMtx);
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
        segAndOrientCh(tempFgMask, tempChFgMask, procFrameArr[ch], impConf.chROIs_[ch],
                       impConf.chWidth_);
        // find fgLocs in center column
        cv::findNonZero(procFrameArr[ch].col(procFrameArr[ch].cols / 2), fgLocs);

        // print fgLocs
        // for (int i = 0; i < fgLocs.size(); ++i)
        //   info("ch: {}, i: {}, fgLocs: {}", ch, i, fgLocs[i]);

        findClusters(fgLocs, fgClstrs[ch]);
        for (auto &fgChClstrs : fgClstrs)
          for (auto &clstrLoc : fgChClstrs)
            clstrLoc *= -1;

        // print each cluster
        // for (int i = 0; i < fgClstrs[ch].size(); ++i)
        //   info("ch: {}, i: {}, fgClstrs: {}", ch, i, fgClstrs[ch][i]);
      }

      fgClstrsFull = findOtherClstrs(fgClstrs);

      // print each cluster
      // for (int ch = 0; ch < impConf.numChs_; ++ch)
      //   for (int i = 0; i < fgClstrsFull[ch].size(); ++i)
      //     info("ch: {}, i: {}, fgClstrsFull: {}", ch, i, fgClstrsFull[ch][i]);

      for (int ch = 0; ch < impConf.numChs_; ++ch) {
        y1[ch] = minDist(fgClstrsFull[ch], yPrev1[ch]);
        y2[ch] = minDist(fgClstrsFull[ch], yPrev2[ch]);
        yPrev1[ch] = y1[ch];
        yPrev2[ch] = y2[ch];
      }

      rstOnZeroCross();

      y.clear();
      y.insert(y.end(), y1.begin(), y1.end());
      y.insert(y.end(), y2.begin(), y2.end());
      procData->push(y);
      procFrameBuf.set(procFrameArr);

      // print y1 and y2
      for (int ch = 0; ch < impConf.numChs_; ++ch)
        info("ch: {}, y1: {}, y2: {}", ch, y1[ch], y2[ch]);

      stopTime = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stopTime - startTime);
      // info("imProc duration: {}", duration.count());
    } catch (cv::Exception &e) {
      error("Message: {}", e.what());
      error("Type: {}", type_name<decltype(e)>());
    }
  }
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
  if (fgLocs.empty()) {
    return;
  }

  int clusterStart = fgLocs[0].y;
  int clusterEnd = fgLocs[0].y;
  for (size_t i = 1; i < fgLocs.size(); ++i) {
    if (fgLocs[i].y == clusterEnd + 1) {
      // This point is part of the current cluster.
      clusterEnd = fgLocs[i].y;
    } else {
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

bool anyZeroCross(const std::vector<double> &vec1, const std::vector<double> &vec2) {
  if (vec1.size() != vec2.size())
    throw std::invalid_argument("Vectors are not the same size");

  auto it = std::mismatch(vec1.begin(), vec1.end(), vec2.begin(),
                          [](double val1, double val2) { return (val1 >= 0) == (val2 >= 0); });

  return it.first != vec1.end();
}

bool anyNonZero(event_bus &eb, std::size_t start, std::size_t end) {
  for (std::size_t i = start; i < end; i++)
    if (eb.r[i] == 0)
      return false;
  return true;
}

void ImProc::rstOnZeroCross() {
  std::size_t rSize = sizeof(currEv.r) / sizeof(currEv.r[0]);

  // if y1 is being controlled and any y2 crosses zero, reset y2
  if (anyNonZero(currEv, 0, rSize / 2))
    if (anyZeroCross(y2, yPrev2))
      std::fill(y2.begin(), y2.end(), 0);

  // if y2 is being controlled and any y1 crosses zero, reset y1
  if (anyNonZero(currEv, rSize / 2, rSize))
    if (anyZeroCross(y1, yPrev1))
      std::fill(y1.begin(), y1.end(), 0);
}

bool ImProc::started() { return startedImProc; }

cv::Mat ImProc::getProcFrame(int idx) { return procFrameBuf.get()[idx]; }

void ImProc::clearProcData() { procData->clear(); }
