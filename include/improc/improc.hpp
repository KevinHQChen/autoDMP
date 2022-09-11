#pragma once

#include "improc/imcap.hpp"
#include "improc/improc.config.hpp"
#include "util/util.hpp"

#define NUM_TEMPLATES 4

struct pose {
  int rotation;
  cv::Point matchLoc;
};

/* Helper functions */
void rotateMat(cv::Mat &src, cv::Mat &dst, double angle);
/* Helper functions */

class ImProc {
  ordered_value conf;
  ordered_value &imProcConf;
  ImCap *imCap = nullptr;
  std::vector<QueueFPS<cv::Mat> *> tempResultQueueArr, procFrameQueueArr;
  cv::Mat preFrame{0, 0, CV_16UC1}, tempFrame{0, 0, CV_16UC1}, tempPreFrame{0, 0, CV_16UC1},
      tempProcFrame{0, 0, CV_16UC1};

  std::vector<int> compParams;

  std::atomic<bool> startedImProc{false};
  std::thread procThread;

  // Called within imProcThread context
  void start();

public:
  // channel/template setup vars
  cv::Rect tmplBBox;
  std::array<cv::Mat, NUM_TEMPLATES> templateImg;
  ImProcConfig impConf;

  ImProc(ImCap *imCap);
  ~ImProc();

  // load/save template images, channel bounding boxes into imProcConfig
  void loadConfig(std::string configPath);
  void saveConfig();

  void startProcThread();
  void stopProcThread();
  bool started();

  std::vector<cv::Mat> getTempFrames();
  cv::Mat getProcFrame(int idx);
  void clearTempFrameQueues();
  void clearProcFrameQueues();
  void setTemplates(std::string tmplSrc);
};
