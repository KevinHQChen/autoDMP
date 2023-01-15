#pragma once

#include "gui/components/button.hpp"
#include "gui/components/dropdown.hpp"
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

  std::string excitationSignalType_;
  std::vector<std::string> excitationSignalTypes_ = {"sine", "square", "triangle", "sawtooth"};
  float excitationSignalMaxValues_[3] = {1.0f, 1.0f, 1.0f};
  float excitationSignalMinValues_[3] = {0.0f, 0.0f, 0.0f};
  float controlSignalSetpoints_[3] = {0.5f, 0.5f, 0.5f};
};

} // namespace gui
