#include "improc/improc.hpp"
#include <numeric>

void findCombs(const std::vector<std::vector<double>> &sets, std::vector<int> &indices,
               int setIndex, std::function<void(const std::vector<int> &)> callback) {
  if (setIndex == sets.size()) {
    callback(indices);
    return;
  }

  for (int i = 0; i < sets[setIndex].size(); ++i) {
    indices[setIndex] = i;
    findCombs(sets, indices, setIndex + 1, callback);
  }
}

void findCombinations(int offset, int j, int m, std::vector<int> &combination,
                      std::vector<std::vector<int>> &combinations) {
  if (j == 0) {
    combinations.push_back(combination);
    return;
  }
  for (int i = offset; i <= m - j; ++i) {
    combination.push_back(i);
    findCombinations(i + 1, j - 1, m, combination, combinations);
    combination.pop_back();
  }
}

double minDist(std::vector<double>& vec, double value) {
    return *std::min_element(vec.begin(), vec.end(), [value](double a, double b) {
        return std::abs(a - value) < std::abs(b - value);
    });
}

std::vector<std::vector<double>> computeNewClstrs(const std::vector<std::vector<double>> &clstrs) {
  int n = clstrs.size();
  std::vector<std::vector<double>> newClstrs(n);

  for (int i = 0; i < n; ++i) {
    std::cout << "##For clstrs " << i << ":\n";
    // store all non-empty elements of `newClstrs` except `newClstrs[i]` in `otherClstrs`
    std::vector<std::vector<double>> otherClstrs(clstrs);
    otherClstrs.erase(otherClstrs.begin() + i);
    otherClstrs.erase(std::remove_if(otherClstrs.begin(), otherClstrs.end(),
                                     [](const std::vector<double> &v) { return v.empty(); }),
                      otherClstrs.end());
    int m = otherClstrs.size();
    std::cout << m << " otherClstrs selected\n";

    for (int j = 1; j < m+1; ++j) {
      std::cout << j << " selected channels\n";
      // find all distinct combinations of j indices from all indices of `otherClstrs`
      std::vector<std::vector<int>> clstrIndices;
      std::vector<int> clstrIdx;
      findCombinations(0, j, m, clstrIdx, clstrIndices);
      int l = clstrIndices.size();
      std::cout << l << " set combinations found\n";
      for (int i = 0; i < l; ++i) {
        std::cout << "clstrIndices[" << i << "] = {";
        for (int j = 0; j < clstrIndices[i].size(); ++j) {
          std::cout << clstrIndices[i][j];
          if (j != clstrIndices[i].size() - 1)
            std::cout << ", ";
        }
        std::cout << "}\n";
      }

      for (int k = 0; k < l; ++k) {
        std::cout << "combination " << k << ", ";
        // given the k-th combination of j indices, store the elements of `otherClstrs` at the
        // corresponding indices in `selClstrs`
        std::vector<std::vector<double>> selClstrs;
        for (int x = 0; x < j; ++x) {
          std::cout << "selClstrs idx " << clstrIndices[k][x];
          selClstrs.push_back(otherClstrs[clstrIndices[k][x]]);
          for (int i = 0; i < otherClstrs[clstrIndices[k][x]].size(); ++i)
            std::cout << ", val " << otherClstrs[clstrIndices[k][x]][i];
          std::cout << "\n";
        }

        // find all distinct combinations of j elements from `selClstrs`
        std::vector<int> indices(j, 0);
        findCombs(selClstrs, indices, 0, [&](const std::vector<int> &indices) {
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

      fgClstrsFull = computeNewClstrs(fgClstrs);

      for (int ch = 0; ch < impConf.numChs_; ++ch) {
        y1[ch] = minDist(fgClstrsFull[ch], yPrev1[ch]);
        y2[ch] = minDist(fgClstrsFull[ch], yPrev2[ch]);
        yPrev1[ch] = y1[ch];
        yPrev2[ch] = y2[ch];
      }

      y.clear();
      y.insert(y.end(), y1.begin(), y1.end());
      y.insert(y.end(), y2.begin(), y2.end());
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
