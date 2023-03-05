#pragma once

#include "gui/components/button.hpp"
#include "gui/components/slider.hpp"
#include "window.hpp"

#include "imcap/imcap.hpp"
#include "improc/improc.hpp"
#include "immvision/immvision.h"

namespace gui {

class ImProcWindow : public Window {
  ImmVision::ImageParams immvisionParams;
  bool imcapVisible_{false}, improcVisible_{false};

  std::shared_ptr<ImCap> imCap_;
  std::shared_ptr<ImProc> imProc_;

  cv::Mat tmpMat, rawMat;

public:
  ImProcWindow(std::shared_ptr<ImCap> imCap, std::shared_ptr<ImProc> imProc);
  ~ImProcWindow();
  void render() override;
};

} // namespace gui
