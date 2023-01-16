#pragma once

#include "ctrl/supervisor.hpp"
#include "gui/components/button.hpp"
#include "gui/components/dropdown.hpp"
#include "gui/components/sliderfloatarray.hpp"
#include "gui/windows/sysid_window.hpp"
#include "window.hpp"

namespace gui {

class SysIdSetupWindow : public Window {
public:
  SysIdSetupWindow(std::shared_ptr<SysIdWindow> sysIDWindow, std::shared_ptr<Supervisor> sv);
  ~SysIdSetupWindow();
  void render() override;

private:
  void previewExcitationSignal();
  void toggleExcitationSignal();
  void sendExcitationSignal();
  void stopExcitationSignal();
  void clearCtrlDataQueue();

  std::shared_ptr<Supervisor> sv_;

  std::weak_ptr<SysIdWindow> sysIdWindow_;
  std::unique_ptr<Button> excitationSignalPreviewBtn_;
  std::unique_ptr<Button> toggleExcitationSignalBtn_;
  std::unique_ptr<Button> sendExcitationSignalBtn_;
  std::unique_ptr<Button> stopExcitationSignalBtn_;
  std::unique_ptr<Button> clearDataBtn_;
  std::unique_ptr<Dropdown> excitationSignalDropdown_;
  std::unique_ptr<SliderFloatArray> minValSlider_;
  std::unique_ptr<SliderFloatArray> maxValSlider_;
  std::unique_ptr<SliderFloatArray> urefSlider_;

  std::vector<std::string> excitationSignalTypes_ = {"sine", "square", "triangle", "sawtooth"};
};

} // namespace gui
