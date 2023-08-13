#pragma once

#include "gui/components/button.hpp"
#include "gui/components/slider.hpp"
#include "window.hpp"

#include "pump/pump.hpp"

namespace gui {

class PumpWindow : public Window {
  void setPumps();
  void syncPumps();
  void setValves();
  void setFreq();
  void setMaxOutput();
  void resetPump();

  Pump *pp_;
  std::string pumpUnits_;

  std::unique_ptr<Slider<int>> maxOutputSlider_;
  std::unique_ptr<SliderArray<float>> outputSlider_;
  std::unique_ptr<Toggle> syncToggle_;
  std::unique_ptr<Slider<int>> freqSlider_;
  std::unique_ptr<ToggleArray> valveToggle_;
  std::unique_ptr<Toggle> controlToggle_;
  std::unique_ptr<Button> resetBtn_;

public:
  PumpWindow(Pump *pp);
  ~PumpWindow();
  void render() override;
};

} // namespace gui
