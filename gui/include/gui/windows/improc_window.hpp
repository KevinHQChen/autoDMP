#pragma once

#include "gui/components/button.hpp"
#include "gui/components/image.hpp"
#include "gui/components/slider.hpp"
#include "window.hpp"

#include "imcap/imcap.hpp"
#include "improc/improc.hpp"

#define NUM_CHANS 3

namespace gui {

class ImProcWindow : public Window {
  bool improcSetupVisible_{false}, improcVisible_{false};

  // ImGuiWindowFlags imCapFlags = 0;
  // ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;

  std::shared_ptr<ImCap> imCap_;
  std::shared_ptr<ImProc> imProc_;

  std::unique_ptr<Toggle> imProcSetupToggle_;
  std::unique_ptr<IMMImage> rawImage_, procImage_, tmplImage1_, tmplImage2_;

  std::array<std::unique_ptr<IMMImage>, NUM_CHANS> chImages_;
  std::array<std::unique_ptr<IMMImage>, NUM_TEMPLATES> tmplImages_;

public:
  ImProcWindow(std::shared_ptr<ImCap> imCap, std::shared_ptr<ImProc> imProc);
  ~ImProcWindow();
  void render() override;
};

} // namespace gui
