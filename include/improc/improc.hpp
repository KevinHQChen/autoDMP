#pragma once

#include "improc/imcap.hpp"
#include "util/util.hpp"

#define NUM_TEMPLATES 4

struct pose {
  int rotation;
  cv::Point matchLoc;
};

struct ChannelPose {
  std::vector<int> rotAngle;
  std::vector<cv::Rect> chanBBox;
  std::vector<cv::Rect> rotChanBBox;
};

class ImProc {
  ordered_value conf;
  ordered_value &imProcConf;
  ImCap *imCap = nullptr;
  std::vector<QueueFPS<cv::Mat> *> tempResultQueueArr, procFrameQueueArr;
  cv::Mat preFrame{0, 0, CV_16UC1};
  std::array<cv::Mat, NUM_TEMPLATES> templateImg;
  ChannelPose chanPose;

  std::vector<int> compParams;

  std::atomic<bool> startedImProc{false}, startedImProcSetup{false};
  std::thread procThread, setupThread;

  // Called within imProcSetupThread context
  void startSetup();

  // Called within imProcThread context
  void start();

public:
  ImProc(ImCap *imCap);
  ~ImProc();

  void startSetupThread();
  void stopSetupThread();
  bool startedSetup();

  void startProcThread();
  void stopProcThread();
  bool started();

  std::vector<cv::Mat> getTempFrames();
  std::vector<cv::Mat> getProcFrames();
};
