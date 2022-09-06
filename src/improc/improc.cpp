#include "improc/improc.hpp"

ImProc::ImProc(ImCap *imCap)
    : conf(toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml")),
      imProcConf(toml::find(conf, "improc")), imCap(imCap),
      tempResultQueueArr({new QueueFPS<cv::Mat>("tempResultsQueue1.txt"),
                          new QueueFPS<cv::Mat>("tempResultsQueue2.txt"),
                          new QueueFPS<cv::Mat>("tempResultsQueue3.txt")}),
      procFrameQueueArr({new QueueFPS<cv::Mat>("procFramesQueue1.txt"),
                         new QueueFPS<cv::Mat>("procFramesQueue2.txt"),
                         new QueueFPS<cv::Mat>("procFramesQueue3.txt")}) {
  // save images with proper format PNG, CV_16UC1
  compParams.push_back(cv::IMWRITE_PNG_COMPRESSION);
  compParams.push_back(0);
}

ImProc::~ImProc() {
  stopSetupThread();
  stopProcThread();
  for (auto &q : procFrameQueueArr)
    delete q;
  for (auto &q : tempResultQueueArr)
    delete q;
}

void ImProc::startSetupThread() {
  info("Starting image processing setup...");
  startedImProcSetup = true;
  imCap->clearPreFrameQueue();
  setupThread = std::thread(&ImProc::startSetup, this);
  setupThread.detach();
}

void ImProc::stopSetupThread() {
  info("Stopping image processing setup...");
  startedImProcSetup = false;
  if (setupThread.joinable())
    setupThread.join();
  imCap->clearPreFrameQueue();
}

void ImProc::startSetup() {
  while(true) {
  if (readFromFile)
    // TODO read channel/template config from file
  if (writeToFile)
    // TODO write channel/template config to file
  }

  if (toml::get<std::string>(imProcConf["tmplSrc"]) == "fromFile") {
  } else {
    // TODO grab template from user-selected frame
  }
}

bool ImProc::startedSetup() { return startedImProcSetup; }

void ImProc::startProcThread() {
  info("Starting image processing...");
  startedImProc = true;
  imCap->clearPreFrameQueue();
  procThread = std::thread(&ImProc::start, this);
  procThread.detach();
}

void ImProc::stopProcThread() {
  info("Stopping image processing...");
  startedImProc = false;
  if (procThread.joinable())
    procThread.join();
  imCap->clearPreFrameQueue();
  for (auto &q : procFrameQueueArr)
    q->clear();
  for (auto &q : tempResultQueueArr)
    q->clear();
}

void ImProc::start() {
  startedImProc = true;
  while (startedImProc) {
    continue;
    // TODO implement image processing
    // for (auto &q : procFrameQueueArr) {
    //   if (q->size() > 0) {
    //     cv::Mat frame = q->pop();
    //     cv::Mat result;
    //     cv::cvtColor(frame, result, cv::COLOR_BGR2GRAY);
    //     tempResultQueueArr[q->getId()]->push(result);
    //   }
    // }
  }
}

bool ImProc::started() { return startedImProc; }

std::vector<cv::Mat> ImProc::getTempFrames() {
  std::vector<cv::Mat> tempFrames;
  for (auto &q : tempResultQueueArr)
    tempFrames.push_back(q->get());
  return tempFrames;
}

std::vector<cv::Mat> ImProc::getProcFrames() {
  std::vector<cv::Mat> procFrames;
  for (auto &q : procFrameQueueArr)
    if (!q->empty())
      procFrames.push_back(q->get());
  return procFrames;
}
