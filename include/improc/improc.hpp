#pragma once

#include "improc/imcap.hpp"
#include "util/util.hpp"

#define NUM_TEMPLATES 4

struct pose {
  int rotation;
  cv::Point matchLoc;
};

class ChannelPose {
  mutable std::mutex mtx;
public:
  std::vector<int> rotAngle;
  cv::Point junction;
  cv::Rect bbox; // x, y, width, height
  int chanWidth;
  std::vector<cv::Rect> chanBBox;
  std::vector<cv::Rect> rotChanBBox;

  ChannelPose() : rotAngle{135, -135, 0}, junction(cv::Point(200, 200)), bbox(cv::Rect(0, 0, junction.x*2, junction.y*2)), chanWidth(50) {
    // generate sensible defaults
    chanBBox.push_back(cv::Rect(junction.x, junction.y, bbox.width/2, bbox.height/2));
    chanBBox.push_back(cv::Rect(0, junction.y, bbox.width/2, bbox.height/2));
    chanBBox.push_back(cv::Rect(junction.x-chanWidth/2, 0, chanWidth, bbox.height/2));

    rotChanBBox.push_back(cv::Rect(bbox.height/4*1.414, 0, chanWidth, bbox.height/2*1.414));
    rotChanBBox.push_back(cv::Rect(bbox.height/4*1.414, 0, chanWidth, bbox.height/2*1.414));
    rotChanBBox.push_back(cv::Rect(0,0,0,0));
  }

  ChannelPose(const ChannelPose& other) {
    rotAngle = other.rotAngle;
    junction = other.junction;
    bbox = other.bbox;
    chanWidth = other.chanWidth;
    chanBBox = other.chanBBox;
    rotChanBBox = other.rotChanBBox;
  }

  ChannelPose &operator=(const ChannelPose &other) {
    // https://stackoverflow.com/a/29609593
    std::lock(mtx, other.mtx);
    std::lock_guard<std::mutex> lhs_lk(mtx, std::adopt_lock);
    std::lock_guard<std::mutex> rhs_lk(other.mtx, std::adopt_lock);

    rotAngle = other.rotAngle;
    junction = other.junction;
    bbox = other.bbox;
    chanWidth = other.chanWidth;
    chanBBox = other.chanBBox;
    rotChanBBox = other.rotChanBBox;
    return *this;
  }

  // ChannelPose &operator=(ChannelPose &&chanPoseInstance) {
  //   std::lock_guard<std::mutex> lockGuard(mtx);
  //   rotAngle = std::move(chanPoseInstance.rotAngle);
  //   chanBBox = std::move(chanPoseInstance.chanBBox);
  //   rotChanBBox = std::move(chanPoseInstance.rotChanBBox);
  //   junction = std::move(chanPoseInstance.junction);
  //   bbox = std::move(chanPoseInstance.bbox);
  //   chanWidth = std::move(chanPoseInstance.chanWidth);
  //   return *this;
  // }
};

/* Helper functions */
void rotateMat(cv::Mat& src, cv::Mat& dst, double angle);
/* Helper functions */

class ImProc {
  ordered_value conf;
  ordered_value &imProcConf;
  ImCap *imCap = nullptr;
  std::vector<QueueFPS<cv::Mat> *> tempResultQueueArr, procFrameQueueArr;
  cv::Mat preFrame{0, 0, CV_16UC1}, tempFrame{0, 0, CV_16UC1}, tempPreFrame{0, 0, CV_16UC1}, tempProcFrame{0, 0, CV_16UC1};

  // channel/template setup vars
  std::array<cv::Mat, NUM_TEMPLATES> templateImg;
  ChannelPose chanPose;
  cv::Mat currChan;

  std::vector<int> compParams;

  std::atomic<bool> startedImProc{false};
  std::thread procThread;

  // Called within imProcThread context
  void start();

public:
  ImProc(ImCap *imCap);
  ~ImProc();

  // load template images into templateImg, channel bounding boxes into chanPose
  void loadConfig(std::string configPath);

  // void startSetupThread();
  // void stopSetupThread();
  // bool startedSetup();

  void startProcThread();
  void stopProcThread();
  bool started();

  std::vector<cv::Mat> getTempFrames();
  cv::Mat getProcFrame(int idx);
  void clearTempFrameQueues();
  void clearProcFrameQueues();
};
