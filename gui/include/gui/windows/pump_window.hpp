#pragma once

#include "gui/components/button.hpp"
#include "gui/components/slider.hpp"
#include "window.hpp"

#include "pump/pump.hpp"
#include <ctrl/state/state.hpp>

#define NUM_PUMPS 4

namespace gui {

class PumpWindow : public Window {
  void setPumps();
  void syncPumps();
  void setValves();
  void setFreq();
  void setMaxVoltage();
  void resetPump();

  std::shared_ptr<Pump> pp_;

  std::unique_ptr<Slider<int>> maxVoltageSlider_;
  std::unique_ptr<SliderArray<int>> voltageSlider_;
  std::unique_ptr<Toggle> syncToggle_;
  std::unique_ptr<Slider<int>> freqSlider_;
  std::unique_ptr<ToggleArray> valveToggle_;
  std::unique_ptr<Toggle> controlToggle_;
  std::unique_ptr<Button> resetBtn_;

public:
  PumpWindow(std::shared_ptr<Pump> pp);
  ~PumpWindow();
  void render() override;
};

} // namespace gui
