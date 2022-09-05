#include "improc/improc.hpp"

ImProc::ImProc(ImCap *imCap)
    : conf(toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml")),
      imProcConf(toml::find(conf, "improc")),
      imCap(imCap), tempResultQueueArr({new QueueFPS<cv::Mat>("tempResultsQueue1.txt"),
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
  stopImProcThread();
  for (auto &q : procFrameQueueArr)
    delete q;
  for (auto &q : tempResultQueueArr)
    delete q;
}

void ImProc::startImProcThread() {
  info("Starting image processing...");
  imProcThread = std::thread(&ImProc::start, this);
  imProcThread.detach();
}

void ImProc::setupTmplMatch() {
  if(toml::get<std::string>(imProcConf["tmplSrc"]) == "fromFile") {
    // TODO load template from file
  } else {
    // TODO grab template from user-selected frame
  }
}

void ImProc::start() {
  setupTmplMatch();
  startImProc = true;
  while (startImProc) {
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

bool ImProc::started() { return startImProc; }

void ImProc::stopImProcThread() {
  info("Stopping image processing...");
  startImProc = false;
  if (imProcThread.joinable())
    imProcThread.join();
  for (auto &q : procFrameQueueArr)
    q->clear();
  for (auto &q : tempResultQueueArr)
    q->clear();
}

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
