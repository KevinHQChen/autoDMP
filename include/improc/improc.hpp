#pragma once

#include "improc/imcap.hpp"
#include "improc/improc.config.hpp"
#include "util/util.hpp"

struct pose {
  int rotation;
  cv::Point matchLoc;
};

/* Helper functions */
void rotateMat(cv::Mat &src, cv::Mat &dst, double angle);
/* Helper functions */

class ImProc {
  ordered_value conf;
  std::string confPath, dataPath;
  int numChans;
  ImCap *imCap = nullptr;
  QueueFPS<cv::Mat> *procFrameQueuePtr;
  std::vector<QueueFPS<cv::Mat> *> tempResultQueueArr, procFrameQueueArr;
  std::vector<QueueFPS<cv::Point> *> procDataQArr;
  std::vector<int> compParams;

  cv::Mat preFrame{0, 0, CV_16UC1}, tempFrame{0, 0, CV_16UC1}, tempPreFrame{0, 0, CV_16UC1},
      tempProcFrame{0, 0, CV_16UC1};

  // for tmpl matching
  cv::Mat tmplFrames[NUM_TEMPLATES];
  cv::Mat tempResultFrame[NUM_TEMPLATES];
  cv::Point minLoc, maxLoc;
  std::optional<cv::Point> currMaxLoc;
  std::vector<cv::Point> maxLocs;
  double minVal, maxVal;

  std::atomic<bool> startedImProc{false}, startedSetup{false};
  std::thread procThread;

  // Called within imProcThread context
  void start();

public:
  ImProcConfig impConf;

  ImProc(ImCap *imCap);
  ~ImProc();

  std::atomic<double> tmplThres{0.75};

  // load/save template images, channel bounding boxes into imProcConfig
  void loadConfig();
  void saveConfig();

  void startProcThread();
  void stopProcThread();
  bool started();
  void setSetupStatus(bool status);

  std::vector<cv::Mat> getTempFrames();
  cv::Mat getProcFrame(int idx);
  cv::Mat getProcFrame();
  void clearTempFrameQueues();
  void clearProcFrameQueues();
};
