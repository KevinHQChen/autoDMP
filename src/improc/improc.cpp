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
  // stopSetupThread();
  stopProcThread();
  for (auto &q : procFrameQueueArr)
    delete q;
  for (auto &q : tempResultQueueArr)
    delete q;
}

// load channel/template images, bounding boxes from file
void ImProc::loadConfig(std::string configPath) {
  // cv::FileStorage chanPoseFile(configPath + "chanPose.yml", cv::FileStorage::READ);
  // int angle;
  // cv::Rect BBox;
  // cv::Rect rotBBox;
  // for (int i = 0; i < toml::get<int>(conf["improc"]["numChans"]); i++) {
  //   if (chanPoseFile.isOpened()) {
  //     chanPoseFile["rotAngle" + std::to_string(i)] >> angle;
  //     chanPose.rotAngle.push_back(angle);
  //     chanPoseFile["chanBBox" + std::to_string(i) + "x"] >> BBox.x;
  //     chanPoseFile["chanBBox" + std::to_string(i) + "y"] >> BBox.y;
  //     chanPoseFile["chanBBox" + std::to_string(i) + "width"] >> BBox.width;
  //     chanPoseFile["chanBBox" + std::to_string(i) + "height"] >> BBox.height;
  //     chanPose.chanBBox.push_back(BBox);
  //     chanPoseFile["rotChanBBox" + std::to_string(i) + "x"] >> rotBBox.x;
  //     chanPoseFile["rotChanBBox" + std::to_string(i) + "y"] >> rotBBox.y;
  //     chanPoseFile["rotChanBBox" + std::to_string(i) + "width"] >> rotBBox.width;
  //     chanPoseFile["rotChanBBox" + std::to_string(i) + "height"] >> rotBBox.height;
  //     chanPose.rotChanBBox.push_back(rotBBox);
  //   }
  //   currChan = firstImage(chanPose.chanBBox[i]).clone(); // save cropped channel as separate image
  //   if (chanPose.rotAngle[i] != 0) {
  //     rotateMat(currChan, currChan, chanPose.rotAngle[i]); // rotate channel to straighten
  //     chans[i] = currChan(chanPose.rotChanBBox[i])
  //                    .clone(); // crop straightened channel from angled channel and use this image
  //                              // to select templates from
  //   } else
  //     chans[i] = currChan;
  // }

  // // load channel bounding boxes
  // for (auto &chan : chanConf) {
  //   std::string chanName = toml::find<std::string>(chan.second, "name");
  //   std::vector<int> chanPoseVec = toml::find<std::vector<int>>(chan.second, "pose");
  //   cv::Rect chanPose(chanPoseVec[0], chanPoseVec[1], chanPoseVec[2], chanPoseVec[3]);
  //   chanPose[chanName] = chanPose;
  // }
}

void ImProc::startProcThread() {
  if (!started()) {
    info("Starting image processing...");
    startedImProc = true;
    imCap->clearPreFrameQueue();
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
  ChannelPose currPose;
  while (started()) {
    currPose = chanPose;
    preFrame = imCap->getPreFrame();
    // if(toml::get<std::string>(conf["cam"]["source"]) == "Webcam") preFrame.convertTo(preFrame, CV_8UC1, 255.0/65535);
    if (!preFrame.empty()) {
      info("currPose bbox: {}", currPose.bbox);
      tempFrame = preFrame(currPose.bbox);
      int idx = 0;
      for (auto &q : procFrameQueueArr) {
        // use currPose to crop preFrame
        info("currPose chanBBox: {}", currPose.chanBBox[idx]);
        tempPreFrame = tempFrame(currPose.chanBBox[idx]);
        if (currPose.rotAngle[idx] != 0) {
          rotateMat(tempPreFrame, tempProcFrame, currPose.rotAngle[idx]);
          tempPreFrame = tempProcFrame(currPose.rotChanBBox[idx]);
        }
        tempProcFrame = tempPreFrame;
        // TODO do template matching here
        // info("tempPreFrame size: {}", tempPreFrame.size());
        // info("tempProcFrame size: {}", tempProcFrame.size());
        // info("currPose rotChanBBox: {}", currPose.rotChanBBox[idx]);
        //
        q->push(tempProcFrame);
        idx++;
      }
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
