#pragma once

#include "gui/components/button.hpp"
// #include "gui/components/image.hpp"
#include "gui/components/slider.hpp"
#include "gui/guiframe.hpp"
#include "window.hpp"

#include "imcap/imcap.hpp"
#include "improc/improc.hpp"

namespace gui {

class ImProcWindow : public Window {
  const int numChans_ = toml::get<int>(Config::conf["improc"]["numChans"]);

  bool improcSetupVisible_{false}, improcVisible_{false};

  // ImGuiWindowFlags imCapFlags = 0;
  // ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;

  ImCap *imCap_;
  ImProc *imProc_;

  std::unique_ptr<Toggle> imProcSetupToggle_;

  // std::unique_ptr<IMMImage> rawImage_, procImage_;
  // std::array<std::unique_ptr<IMMImage>, numChans_> chImages_;
  // std::array<std::unique_ptr<IMMImage>, NUM_TEMPLATES> tmplImages_;

  GUIFrame rawFrame, preFrame;
  GUIFrame procGUIFrames[3];
  std::array<GUIFrame, NUM_TEMPLATES> tmplGUIFrames;
  std::vector<cv::Mat> procFrames;
  std::vector<int> procWidths, procHeights;

public:
  ImProcWindow(ImCap *imCap, ImProc *imProc);
  ~ImProcWindow();
  void render() override;
};

} // namespace gui
