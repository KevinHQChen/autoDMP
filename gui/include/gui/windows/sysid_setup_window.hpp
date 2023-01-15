#pragma once

#include "gui/components/button.hpp"
#include "gui/components/dropdown.hpp"
#include "gui/components/sliderfloatarray.hpp"
#include "window.hpp"

namespace gui {

class SysIdSetupWindow : public Window {
public:
  SysIdSetupWindow();
  ~SysIdSetupWindow();
  void render() override;

private:
  void previewExcitationSignal();
  void toggleExcitationSignal();
  void sendExcitationSignal();

  std::unique_ptr<Button> excitationSignalPreviewButton_;
  std::unique_ptr<Button> toggleExcitationSignalButton_;
  std::unique_ptr<Button> sendExcitationSignalButton_;
  std::unique_ptr<Dropdown> excitationSignalDropdown_;
  std::unique_ptr<SliderFloatArray> minValSlider_;
  std::unique_ptr<SliderFloatArray> maxValSlider_;
  std::unique_ptr<SliderFloatArray> urefSlider_;

  std::vector<std::string> excitationSignalTypes_ = {"sine", "square", "triangle", "sawtooth"};
};

} // namespace gui
