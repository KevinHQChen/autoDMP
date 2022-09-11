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
  stopProcThread();
  for (auto &q : procFrameQueueArr)
    delete q;
  for (auto &q : tempResultQueueArr)
    delete q;
}

// load channel/template images, bounding boxes from file
void ImProc::loadConfig(std::string configPath) {
  ordered_value v =
      toml::parse<toml::discard_comments, tsl::ordered_map>(configPath + "config.toml");
  impConf.from_toml(v);
}

void ImProc::saveConfig() {
  ordered_value newImProcConfig(impConf);
  std::ofstream out(toml::get<std::string>(imProcConf["path"]) + "config.toml");
  out << toml::format(newImProcConfig);
  out.close();
}

void ImProc::startProcThread() {
  if (!started()) {
    info("Starting image processing...");
    startedImProc = true;
    imCap->clearPreFrameQueue();
    if (toml::get<std::string>(imProcConf["confSrc"]) == "fromFile")
      loadConfig(toml::get<std::string>(imProcConf["path"]));
    procThread = std::thread(&ImProc::start, this);
    procThread.detach();
  }
}

void ImProc::stopProcThread() {
  if (started()) {
    info("Stopping image processing...");
    startedImProc = false;
    if (procThread.joinable())
      procThread.join();
    imCap->clearPreFrameQueue();
    this->clearProcFrameQueues();
    this->clearTempFrameQueues();
  }
}

void ImProc::start() {
  while (started()) {
    preFrame = imCap->getPreFrame();
    try {
      if (!preFrame.empty()) {
        tempFrame = preFrame(impConf.getBBox());
        for (unsigned long i = 0; i < procFrameQueueArr.size(); i++) {
          // use chanPose to crop preFrame
          tempPreFrame = tempFrame(impConf.getChanBBox()[i]);
          if (impConf.getRotAngle()[i] != 0) {
            rotateMat(tempPreFrame, tempProcFrame, impConf.getRotAngle()[i]);
            tempPreFrame = tempProcFrame(impConf.getRotChanBBox()[i]);
          }
          tempProcFrame = tempPreFrame;
          // TODO do template matching here

          // info("tempPreFrame size: {}", tempPreFrame.size());
          // info("tempProcFrame size: {}", tempProcFrame.size());
          // info("currPose rotChanBBox: {}", currPose.rotChanBBox[idx]);
          //
          procFrameQueueArr[i]->push(tempProcFrame);
        }
      }
    } catch (cv::Exception &e) {
      error("Message: {}", e.what());
      error("Type: {}", type_name<decltype(e)>());
    }
  }
}

bool ImProc::started() { return startedImProc; }

std::vector<cv::Mat> ImProc::getTempFrames() {
  std::vector<cv::Mat> tempFrames;
  for (auto &q : tempResultQueueArr)
    tempFrames.push_back(q->get());
  return tempFrames;
}

cv::Mat ImProc::getProcFrame(int idx) {
  if (!procFrameQueueArr[idx]->empty())
    return procFrameQueueArr[idx]->get();
  return cv::Mat();
}

void ImProc::clearTempFrameQueues() {
  for (auto &q : tempResultQueueArr)
    q->clear();
}

void ImProc::clearProcFrameQueues() {
  for (auto &q : procFrameQueueArr)
    q->clear();
}

void ImProc::setTemplates(std::string tmplSrc) {}
