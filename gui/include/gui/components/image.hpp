#pragma once

#include "imgui.h"
#include <cassert>
#include <functional>
#include <string>

#include "immvision/immvision.h"

namespace gui {

class IMMImage {
  std::string label_;
  float scale_;

public:
  cv::Mat image_;
  ImmVision::ImageParams imParams;

  IMMImage(std::string label, float scale, bool interactive = false, std::string zoomKey = "z")
      : label_(label), scale_(scale) {
    if (!interactive)
      imParams = ImmVision::FactorImageParamsDisplayOnly();
    else {
      imParams = ImmVision::ImageParams();
      imParams.ImageDisplaySize = cv::Size(0, 100); // default size
      imParams.ZoomKey = zoomKey;
    }
    imParams.RefreshImage = true; // this enables live video
  }

  void render(const cv::Mat &image) {
    if (!image.empty())
      image_ = image;

    if (!image_.empty()) {
      imParams.ImageDisplaySize = cv::Size(0, image_.rows * scale_); // scale height of any image
      ImmVision::Image(label_, image_, &imParams);
    }
  }
};

} // namespace gui
