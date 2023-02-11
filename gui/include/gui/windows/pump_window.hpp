#pragma once

#include "gui/components/button.hpp"
#include "gui/components/checkboxarray.hpp"
#include "gui/components/dropdown.hpp"
#include "gui/components/slider.hpp"
#include "pump/pump.hpp"
#include "window.hpp"
#include <ctrl/state/state.hpp>

#define NUM_PUMPS 4

namespace gui {

class PumpWindow : public Window {
  bool controlFlag_{false};

  void setPumps();
  void setValves();
  void setFreq();
  void resetPump();

  std::shared_ptr<Pump> pp_;

  std::unique_ptr<Slider<int>> voltageSlider_;
  std::unique_ptr<Slider<int>> freqSlider_;
  std::unique_ptr<Button> valveOnOff_;
  std::unique_ptr<Button> controlOnOff_;
  std::unique_ptr<Button> resetBtn_;

public:
  PumpWindow(std::shared_ptr<Pump> pp);
  ~PumpWindow();
  void render() override;
};

} // namespace gui
