#pragma once
#include "util/util.hpp"

#define NUM_TEMPLATES 2

// macros to generate from_toml/into_toml functions for cv::Point and cv::Rect
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(cv::Point, x, y)
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(cv::Rect, x, y, width, height)

class ImProcConfig {
public:
  mutable std::mutex rotAngleMtx;
  std::vector<int> rotAngle_;
  mutable std::mutex junctionMtx;
  cv::Point junction_;
  mutable std::mutex bboxMtx;
  cv::Rect bbox_; // x, y, width, height
  mutable std::mutex tmplBBoxMtx;
  cv::Rect tmplBBox_;
  mutable std::mutex tmplMtx;
  std::array<cv::Mat, NUM_TEMPLATES> tmplImg_;
  mutable std::mutex chanWidthMtx;
  int chanWidth_;
  mutable std::mutex chanBBoxMtx;
  std::vector<cv::Rect> chanBBox_;
  mutable std::mutex rotChanBBoxMtx;
  std::vector<cv::Rect> rotChanBBox_;

  ImProcConfig()
      : rotAngle_{135, -135, 0}, junction_(cv::Point(200, 200)),
        bbox_(cv::Rect(0, 0, junction_.x * 3, junction_.y * 2)), tmplBBox_(cv::Rect(5, 5, 90, 90)),
        chanWidth_(100) {
    // generate sensible defaults
    for (int i = 0; i < NUM_TEMPLATES; ++i)
      tmplImg_[i] = cv::Mat::zeros(tmplBBox_.size(), CV_8UC1);
    chanBBox_.push_back(cv::Rect(junction_.x, junction_.y, bbox_.width / 2, bbox_.height / 2));
    chanBBox_.push_back(cv::Rect(0, junction_.y, bbox_.width / 2, bbox_.height / 2));
    chanBBox_.push_back(cv::Rect(junction_.x - chanWidth_ / 2, 0, chanWidth_, bbox_.height / 2));

    rotChanBBox_.push_back(
        cv::Rect(bbox_.height / 4.0 * 1.414, 0, chanWidth_, bbox_.height / 2.0 * 1.414));
    rotChanBBox_.push_back(
        cv::Rect(bbox_.height / 4.0 * 1.414, 0, chanWidth_, bbox_.height / 2.0 * 1.414));
    rotChanBBox_.push_back(cv::Rect(0, 0, 0, 0));
  }

  ImProcConfig(const ImProcConfig &other) {
    rotAngle_ = other.getRotAngle();
    junction_ = other.getJunction();
    bbox_ = other.getBBox();
    tmplBBox_ = other.getTmplBBox();
    tmplImg_ = other.getTmplImg();
    chanWidth_ = other.getChanWidth();
    chanBBox_ = other.getChanBBox();
    rotChanBBox_ = other.getRotChanBBox();
  }

  ImProcConfig &operator=(const ImProcConfig &other) {
    rotAngle_ = other.getRotAngle();
    junction_ = other.getJunction();
    bbox_ = other.getBBox();
    tmplBBox_ = other.getTmplBBox();
    tmplImg_ = other.getTmplImg();
    chanWidth_ = other.getChanWidth();
    chanBBox_ = other.getChanBBox();
    rotChanBBox_ = other.getRotChanBBox();

    // // https://stackoverflow.com/a/29609593
    // std::lock(rotAngleMtx, other.rotAngleMtx);
    // std::lock_guard<std::mutex> lhs_lk1(rotAngleMtx, std::adopt_lock);
    // std::lock_guard<std::mutex> rhs_lk1(other.rotAngleMtx, std::adopt_lock);
    // std::lock(junctionMtx, other.junctionMtx);
    // std::lock_guard<std::mutex> lhs_lk2(junctionMtx, std::adopt_lock);
    // std::lock_guard<std::mutex> rhs_lk2(other.junctionMtx, std::adopt_lock);
    // std::lock(bboxMtx, other.bboxMtx);
    // std::lock_guard<std::mutex> lhs_lk3(bboxMtx, std::adopt_lock);
    // std::lock_guard<std::mutex> rhs_lk3(other.bboxMtx, std::adopt_lock);
    // std::lock(chanWidthMtx, other.chanWidthMtx);
    // std::lock_guard<std::mutex> lhs_lk4(chanWidthMtx, std::adopt_lock);
    // std::lock_guard<std::mutex> rhs_lk4(other.chanWidthMtx, std::adopt_lock);
    // std::lock(chanBBoxMtx, other.chanBBoxMtx);
    // std::lock_guard<std::mutex> lhs_lk5(chanBBoxMtx, std::adopt_lock);
    // std::lock_guard<std::mutex> rhs_lk5(other.chanBBoxMtx, std::adopt_lock);
    // std::lock(rotChanBBoxMtx, other.rotChanBBoxMtx);
    // std::lock_guard<std::mutex> lhs_lk6(rotChanBBoxMtx, std::adopt_lock);
    // std::lock_guard<std::mutex> rhs_lk6(other.rotChanBBoxMtx, std::adopt_lock);
    // rotAngle_ = other.rotAngle_;
    // junction_ = other.junction_;
    // bbox_ = other.bbox_;
    // chanWidth_ = other.chanWidth_;
    // chanBBox_ = other.chanBBox_;
    // rotChanBBox_ = other.rotChanBBox_;
    return *this;
  }

  // enables stuff like:
  // impConf = toml::find<ImProcConfig>(v, "ImProcConfig");
  void from_toml(const ordered_value &v) {
    setRotAngle(toml::find<std::vector<int>>(v, "rotAngle"));
    setJunction(toml::find<cv::Point>(v, "junction"));
    setBBox(toml::find<cv::Rect>(v, "bbox"));
    setTmplBBox(toml::find<cv::Rect>(v, "tmplBBox"));
    setChanWidth(toml::find<int>(v, "chanWidth"));
    setChanBBox(toml::find<std::vector<cv::Rect>>(v, "chanBBox"));
    setRotChanBBox(toml::find<std::vector<cv::Rect>>(v, "rotChanBBox"));
  }

  // enables stuff like:
  // ordered_value v(imProcConfig);
  ordered_value into_toml() const {
    ordered_value v;
    v["rotAngle"] = getRotAngle();
    v["junction"] = getJunction();
    v["bbox"] = getBBox();
    v["tmplBBox"] = getTmplBBox();
    v["chanWidth"] = getChanWidth();
    v["chanBBox"] = getChanBBox();
    v["rotChanBBox"] = getRotChanBBox();
    info(v);
    return v;
  }

  void setRotAngle(const std::vector<int> &rotAngle) {
    std::lock_guard<std::mutex> lock(rotAngleMtx);
    rotAngle_ = rotAngle;
  }

  void setJunction(const cv::Point &junction) {
    std::lock_guard<std::mutex> lock(junctionMtx);
    junction_ = junction;
  }

  void setBBox(const cv::Rect &bbox) {
    std::lock_guard<std::mutex> lock(bboxMtx);
    bbox_ = bbox;
  }

  void setTmplBBox(const cv::Rect &tmplBBox) {
    std::lock_guard<std::mutex> lock(tmplBBoxMtx);
    tmplBBox_ = tmplBBox;
  }

  void setTmplImg(const std::array<cv::Mat, NUM_TEMPLATES> &tmplImg) {
    std::lock_guard<std::mutex> lock(tmplMtx);
    tmplImg_ = tmplImg;
  }

  void setChanWidth(int chanWidth) {
    std::lock_guard<std::mutex> lock(chanWidthMtx);
    chanWidth_ = chanWidth;
  }

  void setChanBBox(const std::vector<cv::Rect> &chanBBox) {
    std::lock_guard<std::mutex> lock(chanBBoxMtx);
    chanBBox_ = chanBBox;
  }

  void setRotChanBBox(const std::vector<cv::Rect> &rotChanBBox) {
    std::lock_guard<std::mutex> lock(rotChanBBoxMtx);
    rotChanBBox_ = rotChanBBox;
  }

  std::vector<int> getRotAngle() const {
    std::lock_guard<std::mutex> lock(rotAngleMtx);
    return rotAngle_;
  }

  cv::Point getJunction() const {
    std::lock_guard<std::mutex> lock(junctionMtx);
    return junction_;
  }

  cv::Rect getBBox() const {
    std::lock_guard<std::mutex> lock(bboxMtx);
    return bbox_;
  }

  cv::Rect getTmplBBox() const {
    std::lock_guard<std::mutex> lock(tmplBBoxMtx);
    return tmplBBox_;
  }

  std::array<cv::Mat, NUM_TEMPLATES> getTmplImg() const {
    std::lock_guard<std::mutex> lock(tmplMtx);
    return tmplImg_;
  }

  int getChanWidth() const {
    std::lock_guard<std::mutex> lock(chanWidthMtx);
    return chanWidth_;
  }

  std::vector<cv::Rect> getChanBBox() const {
    std::lock_guard<std::mutex> lock(chanBBoxMtx);
    return chanBBox_;
  }

  std::vector<cv::Rect> getRotChanBBox() const {
    std::lock_guard<std::mutex> lock(rotChanBBoxMtx);
    return rotChanBBox_;
  }
};

// ImProcConfig(const ImProcConfig &other) {
//   rotAngle = other.rotAngle;
//   junction = other.junction;
//   bbox = other.bbox;
//   chanWidth = other.chanWidth;
//   chanBBox = other.chanBBox;
//   rotChanBBox = other.rotChanBBox;
// }

// ImProcConfig &operator=(const ImProcConfig &other) {
//   // https://stackoverflow.com/a/29609593
//   std::lock(mtx, other.mtx);
//   std::lock_guard<std::mutex> lhs_lk(mtx, std::adopt_lock);
//   std::lock_guard<std::mutex> rhs_lk(other.mtx, std::adopt_lock);

//   rotAngle = other.rotAngle;
//   junction = other.junction;
//   bbox = other.bbox;
//   chanWidth = other.chanWidth;
//   chanBBox = other.chanBBox;
//   rotChanBBox = other.rotChanBBox;
//   return *this;
// }

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
