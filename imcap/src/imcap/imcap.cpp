#include "imcap/imcap.hpp"

ImCap::ImCap()
    : conf(Config::conf), dataPath(toml::get<std::string>(conf["postproc"]["rawDataPath"])),
      cam(new Cam(0, conf)) {
  info("Initializing image capture...");
}

ImCap::~ImCap() {
  info("Terminating image capture...");
  stopThread();
  delete cam;
}

void ImCap::startThread() {
  if (!started()) {
    info("Starting image capture...");
    startedImCap = true;
    cam->start((int)(100 / 1000)); // timerInterval of 100ms
    captureThread = std::thread(&ImCap::start, this);
    captureThread.detach();
  }
}

void ImCap::stopThread() {
  if (started()) {
    info("Stopping image capture...");
    startedImCap = false;
    std::lock_guard<std::mutex> guard(imcapMtx); // wait for thread to finish
    cam->stop();
  }
}

void ImCap::start() {
  while (startedImCap) {
    std::lock_guard<std::mutex> guard(imcapMtx);
    // Timer t("ImCap");
    imCapSuccess = cam->process(currImg);
    if (imCapSuccess) {
      if (toml::get<std::string>(conf["cam"]["source"]) == "Andor")
        currImg.convertTo(currImg, CV_8UC1, 255.0 / 65535);
      rawFrameBuf.set(currImg);
    } else
      continue; // error("cannot read image");
  }
}

bool ImCap::started() { return startedImCap; }

cv::Mat ImCap::getFrame() { return rawFrameBuf.get(); }
