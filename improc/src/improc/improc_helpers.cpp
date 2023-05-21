#include "improc/improc.hpp"

// from https://stackoverflow.com/a/24352524
void rotateMat(cv::Mat &src, cv::Mat &dst, double angle) {
  // get rotation matrix for rotating the image around its center in pixel coordinates
  cv::Point2f center((src.cols - 1) / 2.0, (src.rows - 1) / 2.0);
  cv::Mat rot = cv::getRotationMatrix2D(center, angle, 1.0);
  // determine bounding rectangle, center not relevant
  cv::Rect2f bbox = cv::RotatedRect(cv::Point2f(), src.size(), angle).boundingRect2f();
  // adjust transformation matrix
  rot.at<double>(0, 2) += bbox.width / 2.0 - src.cols / 2.0;
  rot.at<double>(1, 2) += bbox.height / 2.0 - src.rows / 2.0;

  cv::warpAffine(src, dst, rot, bbox.size());
}

// macros to generate from_toml/into_toml functions for cv::Point and cv::Rect
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(cv::Point, x, y)
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(cv::Rect, x, y, width, height)
TOML11_DEFINE_CONVERSION_NON_INTRUSIVE(RotRect, x, y, width, height, chHeight, angle)
ImProcConfig::ImProcConfig() {
  chROIs_.push_back(RotRect(0, 0, 100, 100, 100));
  chROIs_.push_back(RotRect(0, 0, 100, 100, 100));
  chROIs_.push_back(RotRect(0, 0, 100, 100, 100));
  chWidth_ = 100;
  numChs_ = 3;
  bgSubHistory_ = 500;
  bgSubThres_ = 16;
}

ImProcConfig::ImProcConfig(const ImProcConfig &other) {
  chROIs_ = other.getChROIs();
  chWidth_ = other.getChWidth();
  numChs_ = other.getNumChs();
  bgSubHistory_ = other.getBgSubHistory();
  bgSubThres_ = other.getBgSubThres();
}

// enables stuff like:
// impConf = toml::find<ImProcConfig>(v, "ImProcConfig");
void ImProcConfig::from_toml(const ordered_value &v) {
  setChROIs(toml::find<std::vector<RotRect>>(v, "chROIs"));
  chWidth_ = toml::find<int>(v, "chWidth");
  numChs_ = toml::find<int>(v, "numChs");
  bgSubHistory_ = toml::find<int>(v, "bgSubHistory");
  bgSubThres_ = toml::find<double>(v, "bgSubThres");
}

// enables stuff like:
// ordered_value v(imProcConfig);
ordered_value ImProcConfig::into_toml() const {
  ordered_value v;
  v["chROIs"] = getChROIs();
  v["chWidth"] = getChWidth();
  v["numChs"] = getNumChs();
  v["bgSubHistory"] = getBgSubHistory();
  v["bgSubThres"] = getBgSubThres();
  info(v);
  return v;
}

std::vector<RotRect> ImProcConfig::getChROIs() const {
  std::lock_guard<std::mutex> lock(chanROIMtx);
  return chROIs_;
}
int ImProcConfig::getChWidth() const {
  std::lock_guard<std::mutex> lock(chWidthMtx);
  return chWidth_;
}
int ImProcConfig::getNumChs() const {
  std::lock_guard<std::mutex> lock(numChsMtx);
  return numChs_;
}
int ImProcConfig::getBgSubHistory() const {
  std::lock_guard<std::mutex> lock(bgSubHistoryMtx);
  return bgSubHistory_;
}
double ImProcConfig::getBgSubThres() const {
  std::lock_guard<std::mutex> lock(bgSubThresMtx);
  return bgSubThres_;
}

void ImProcConfig::setChROIs(const std::vector<RotRect> &chROIs) {
  std::lock_guard<std::mutex> lock(chanROIMtx);
  chROIs_ = chROIs;
}
void ImProcConfig::setChWidth(int chWidth) {
  std::lock_guard<std::mutex> lock(chWidthMtx);
  chWidth_ = chWidth;
}
void ImProcConfig::setNumChs(int numChs) {
  std::lock_guard<std::mutex> lock(numChsMtx);
  numChs_ = numChs;
}
void ImProcConfig::setBgSubHistory(int bgSubHistory) {
  std::lock_guard<std::mutex> lock(bgSubHistoryMtx);
  bgSubHistory_ = bgSubHistory;
}
void ImProcConfig::setBgSubThres(double bgSubThres) {
  std::lock_guard<std::mutex> lock(bgSubThresMtx);
  bgSubThres_ = bgSubThres;
}
