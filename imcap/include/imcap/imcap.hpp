#pragma once

#include "cam/cam.hpp"
#include "util/util.hpp"

class ImCap {
  ordered_value conf;
  std::string dataPath;
  Cam *cam;
  QueueFPS<cv::Mat> *rawFrameQueue, *preFrameQueue;
  cv::Mat currImg{0, 0, CV_16UC1};
  bool imCapSuccess;
  std::atomic<bool> startedImCap{false};
  std::thread captureThread;

  // Called within captureThread context
  void start();

public:
  ImCap();
  ~ImCap();
  void startCaptureThread();
  void stopCaptureThread();
  cv::Mat getRawFrame();
  cv::Mat getPreFrame();
  void clearRawFrameQueue();
  void clearPreFrameQueue();
  bool started();
};
