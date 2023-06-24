#include "imcap/imcap.hpp"

ImCap::ImCap()
    : conf(Config::conf), dataPath(toml::get<std::string>(conf["postproc"]["rawDataPath"])),
      cam(new Cam(0, conf)), rawFrameBuf(new SharedBuffer<cv::Mat>()) {}

ImCap::~ImCap() {
  stopThread();
  delete rawFrameBuf;
  delete cam;
}

void ImCap::startThread() {
  if (!started()) {
    info("Starting image capture...");
    startedImCap = true;
    captureThread = std::thread(&ImCap::start, this);
    captureThread.detach();
  }
}

void ImCap::stopThread() {
  if (started()) {
    info("Stopping image capture...");
    startedImCap = false;
    if (captureThread.joinable())
      captureThread.join();
    cam->stop();
    rawFrameBuf->clear();
  }
}

void ImCap::start() {
  cam->start((int)(100 / 1000)); // timerInterval of 100ms
  while (startedImCap) {
    // Timer t("ImCap");
    imCapSuccess = cam->process(currImg);
    if (imCapSuccess) {
      if (toml::get<std::string>(conf["cam"]["source"]) == "Andor")
        currImg.convertTo(currImg, CV_8UC1, 255.0 / 65535);
      rawFrameBuf->set(currImg.clone());
    } else
      continue; // error("cannot read image");
  }
}

bool ImCap::started() { return startedImCap; }

cv::Mat ImCap::getFrame() { return rawFrameBuf->get(); }
