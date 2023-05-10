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
  numTmpls_ = 2;
  numChs_ = 3;
  tmplThres_ = 0.8;
}

ImProcConfig::ImProcConfig(const ImProcConfig &other) {
  chROIs_ = other.getChROIs();
  chWidth_ = other.getChWidth();
  numTmpls_ = other.getNumTmpls();
  numChs_ = other.getNumChs();
  tmplThres_ = other.getTmplThres();
}

// enables stuff like:
// impConf = toml::find<ImProcConfig>(v, "ImProcConfig");
void ImProcConfig::from_toml(const ordered_value &v) {
  setChROIs(toml::find<std::vector<RotRect>>(v, "chROIs"));
  chWidth_ = toml::find<int>(v, "chWidth");
  numTmpls_ = toml::find<int>(v, "numTmpls");
  numChs_ = toml::find<int>(v, "numChs");
  tmplThres_ = toml::find<double>(v, "tmplThres");
}

// enables stuff like:
// ordered_value v(imProcConfig);
ordered_value ImProcConfig::into_toml() const {
  ordered_value v;
  v["chROIs"] = getChROIs();
  v["chWidth"] = getChWidth();
  v["numTmpls"] = getNumTmpls();
  v["numChs"] = getNumChs();
  v["tmplThres"] = getTmplThres();
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
int ImProcConfig::getNumTmpls() const {
  std::lock_guard<std::mutex> lock(numTmplsMtx);
  return numTmpls_;
}
int ImProcConfig::getNumChs() const {
  std::lock_guard<std::mutex> lock(numChsMtx);
  return numChs_;
}
double ImProcConfig::getTmplThres() const {
  std::lock_guard<std::mutex> lock(tmplThresMtx);
  return tmplThres_;
}
std::vector<cv::Mat> ImProcConfig::getTmplImg() const {
  std::lock_guard<std::mutex> lock(tmplMtx);
  return tmplImg_;
}

void ImProcConfig::setChROIs(const std::vector<RotRect> &chROIs) {
  std::lock_guard<std::mutex> lock(chanROIMtx);
  chROIs_ = chROIs;
}
void ImProcConfig::setChWidth(int chWidth) {
  std::lock_guard<std::mutex> lock(chWidthMtx);
  chWidth_ = chWidth;
}
void ImProcConfig::setNumTmpls(int numTmpls) {
  std::lock_guard<std::mutex> lock(numTmplsMtx);
  numTmpls_ = numTmpls;
}
void ImProcConfig::setNumChs(int numChs) {
  std::lock_guard<std::mutex> lock(numChsMtx);
  numChs_ = numChs;
}
void ImProcConfig::setTmplThres(double tmplThres) {
  std::lock_guard<std::mutex> lock(tmplThresMtx);
  tmplThres_ = tmplThres;
}
void ImProcConfig::setTmplImg(cv::Mat tmplImg) {
  std::lock_guard<std::mutex> lock(tmplMtx);
  tmplImg_.push_back(tmplImg);
}
void ImProcConfig::setTmplImg(const std::vector<cv::Mat> &tmplImg) {
  std::lock_guard<std::mutex> lock(tmplMtx);
  tmplImg_ = tmplImg;
}
void ImProcConfig::clearTmplImgs() {
  std::lock_guard<std::mutex> lock(tmplMtx);
  tmplImg_.clear();
}
