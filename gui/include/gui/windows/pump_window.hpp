#pragma once

#include "pump/pump.hpp"
#include "gui/components/button.hpp"
#include "gui/components/checkboxarray.hpp"
#include "gui/components/dropdown.hpp"
#include "gui/components/sliderint.hpp"
#include "window.hpp"
#include <ctrl/state/state.hpp>

#define NUM_PUMPS 4

namespace gui {

class PumpWindow : public Window {
  void setPumps();
  void setValves();
  void setFreq();

  std::shared_ptr<Pump> pp_;

  std::unique_ptr<SliderInt> voltageSlider_;
  std::unique_ptr<SliderInt> freqSlider_;
  std::unique_ptr<Button> valveOnOff_;

public:
  PumpWindow(std::shared_ptr<Pump> pp);
  ~PumpWindow();
  void render() override;
};

} // namespace gui
