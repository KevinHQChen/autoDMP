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
