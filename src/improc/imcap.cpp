#include "improc/imcap.hpp"

ImCap::ImCap()
    : conf(toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml")),
      dataPath(toml::get<std::string>(conf["postproc"]["rawDataPath"])), cam(new Cam(0, conf)),
      rawFrameQueuePtr(new QueueFPS<cv::Mat>(dataPath + "rawFramesQueue.txt")),
      preFrameQueuePtr(new QueueFPS<cv::Mat>(dataPath + "preFramesQueue.txt")) {}

ImCap::~ImCap() {
  stopCaptureThread();
  delete preFrameQueuePtr;
  delete rawFrameQueuePtr;
  delete cam;
}

void ImCap::startCaptureThread() {
  if (!started()) {
    info("Starting image capture...");
    startedImCap = true;
    captureThread = std::thread(&ImCap::start, this);
    captureThread.detach();
  }
}

void ImCap::stopCaptureThread() {
  if (started()) {
    info("Stopping image capture...");
    startedImCap = false;
    if (captureThread.joinable())
      captureThread.join();
    cam->stop();
    rawFrameQueuePtr->clear();
    preFrameQueuePtr->clear();
  }
}

void ImCap::start() {
  // start camera (allocate circular buffer to store frames)
  // timerInterval sets the size of buffers needed (min size of 1) multiplied by the framerate
  // (pretty sure) if timerInterval is less than 1000ms, only 1 buffer (i.e .40 frames) is needed
  cam->start((int)(100 / 1000)); // timerInterval of 100ms
  while (startedImCap) {
    // auto startTime = high_resolution_clock::now();
    imCapSuccess = cam->process(currImg);
    if (imCapSuccess) {
      rawFrameQueuePtr->push(currImg.clone());
      if (toml::get<std::string>(conf["cam"]["source"]) == "Andor")
        currImg.convertTo(currImg, CV_8UC1, 255.0 / 65535);
      preFrameQueuePtr->push(currImg.clone());
    } else
      continue; // error("cannot read image");
    // auto stopTime = high_resolution_clock::now();
    // auto duration = duration_cast<milliseconds>(stopTime - startTime);
    // info("imCap duration: {}", duration.count());
  }
}

bool ImCap::started() { return startedImCap; }

cv::Mat ImCap::getRawFrame() {
  if (!rawFrameQueuePtr->empty())
    return rawFrameQueuePtr->get();
  return cv::Mat();
}

cv::Mat ImCap::getPreFrame() {
  if (!preFrameQueuePtr->empty())
    return preFrameQueuePtr->get();
  return cv::Mat();
}

void ImCap::clearRawFrameQueue() { rawFrameQueuePtr->clear(); }

void ImCap::clearPreFrameQueue() { preFrameQueuePtr->clear(); }
