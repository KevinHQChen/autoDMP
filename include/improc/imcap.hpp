#pragma once

#include "cam/cam.hpp"
#include "util/util.hpp"

class ImCap {
  ordered_value conf;
  Cam *cam;
  QueueFPS<cv::Mat> *rawFrameQueuePtr, *preFrameQueuePtr;
  cv::Mat currImg{0, 0, CV_16UC1};
  bool imCapSuccess;
  std::atomic<bool> startImCap{false};
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
  bool started();
};
