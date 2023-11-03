#pragma once

#include "cam/cam.hpp"
#include "util/util.hpp"

class ImCap {
public:
  ImCap(Cam *cam, std::shared_ptr<logger> log);
  ~ImCap();
  void startThread();
  void stopThread();
  cv::Mat getRawFrame();
  cv::Mat getFrame();
  void clearRawFrameQueue();
  void clearPreFrameQueue();
  bool started();
  void saveRawFrameToFile();

private:
  std::shared_ptr<logger> lg;
  ordered_value conf;
  std::string dataPath;
  std::vector<int> compParams;
  Cam *cam_;

  std::mutex imcapMtx;
  std::atomic<bool> startedImCap{false};
  std::thread captureThread;
  SharedBuffer<cv::Mat> rawFrameBuf;

  cv::Mat currImg{0, 0, CV_16UC1};
  bool imCapSuccess;

  /*
  ** start camera (allocate circular buffer to store frames)
  ** timerInterval sets the size of buffers needed (min size of 1) multiplied by the framerate
  ** (pretty sure) if timerInterval is less than 1000ms, only 1 buffer (i.e .40 frames) is needed
  ** (Called within captureThread context)
  */
  void start();
};
