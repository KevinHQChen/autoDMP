#include "imcap/imcap.hpp"

ImCap::ImCap(Cam *cam, std::shared_ptr<logger> log)
    : conf(Config::conf), dataPath(toml::get<std::string>(conf["postproc"]["rawDataPath"])),
      cam_(cam), lg(log) {
  lg->info("Initializing image capture...");

  // save images with proper format PNG, CV_16UC1
  compParams.push_back(cv::IMWRITE_PNG_COMPRESSION);
  compParams.push_back(0);
}

ImCap::~ImCap() {
  lg->info("Terminating image capture...");
  stopThread();
  delete cam_;
}

void ImCap::startThread() {
  if (!started()) {
    lg->info("Starting image capture...");
    startedImCap = true;
    cam_->start((int)(100 / 1000)); // timerInterval of 100ms
    captureThread = std::thread(&ImCap::start, this);
    captureThread.detach();
  }
}

void ImCap::stopThread() {
  if (started()) {
    lg->info("Stopping image capture...");
    startedImCap = false;
    std::lock_guard<std::mutex> guard(imcapMtx); // wait for thread to finish
    cam_->stop();
  }
}

void ImCap::start() {
  while (startedImCap) {
    std::lock_guard<std::mutex> guard(imcapMtx);
    // Timer t("ImCap");
    imCapSuccess = cam_->process(currImg);
    if (imCapSuccess) {
      if (toml::get<std::string>(conf["cam"]["source"]) == "Andor")
        currImg.convertTo(currImg, CV_8UC1, 255.0 / 65535);
      rawFrameBuf.set(currImg);
    } else
      continue; // lg->error("cannot read image");
  }
}

bool ImCap::started() { return startedImCap; }

cv::Mat ImCap::getFrame() { return rawFrameBuf.get(); }

void ImCap::saveRawFrameToFile() {
  try {
    cv::imwrite("output.png", getFrame(), compParams);
  } catch (cv::Exception &e) {
    lg->error("Message: {}", e.what());
    lg->error("Type: {}", type_name<decltype(e)>());
  }
}
