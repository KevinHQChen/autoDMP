#pragma once

#include "imcap/imcap.hpp"
#include "improc/improc.config.hpp"
#include "util/util.hpp"

struct Pose {
  int rot;
  cv::Point loc;
};

/* Helper functions */
void rotateMat(cv::Mat &src, cv::Mat &dst, double angle);
/* Helper functions */

class ImProc {
  ordered_value conf;
  std::string confPath, dataPath;
  int numChans;
  ImCap *imCap;
  QueueFPS<cv::Mat> *procFrameQueuePtr;
  std::vector<QueueFPS<cv::Mat> *> tempResultQueueArr, procFrameQueueArr;
  std::vector<int> compParams;

  cv::Mat preFrame{0, 0, CV_16UC1}, tempFrame{0, 0, CV_16UC1}, tempPreFrame{0, 0, CV_16UC1},
      tempProcFrame{0, 0, CV_16UC1};

  // for tmpl matching
  std::vector<cv::Mat> tmplFrames;
  std::vector<cv::Mat> tempResultFrame;
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
  std::atomic<double> tmplThres;
  std::vector<QueueFPS<Pose> *> procDataQArr;

  ImProc(ImCap *imCap);
  ~ImProc();

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
  cv::Point getProcData(int idx);
  void clearTempFrameQueues();
  void clearProcFrameQueues();
  void clearProcDataQueues();
};
