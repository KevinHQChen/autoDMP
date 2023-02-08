#include "gui/windows/pump_window.hpp"

namespace gui {

PumpWindow::PumpWindow(std::shared_ptr<Pump> pp) : pp_(pp) {
  voltageSlider_ = std::make_unique<SliderInt>("Voltage", 0, 250, "%d V", NUM_PUMPS, ImVec2(36, 200), false, [this]() { setPumps(); });
  freqSlider_ = std::make_unique<SliderInt>("Frequency", 0, 800, "%d Hz\nFreq", 1, ImVec2(36, 200), true, [this]() { setFreq(); });
  valveOnOff_ = std::make_unique<Button>(std::vector<std::string>(NUM_PUMPS, "Open\nValve"), [this]() { setValves(); }, NUM_PUMPS, ImVec2(36, 0));
}

PumpWindow::~PumpWindow() { voltageSlider_.reset(); }

void PumpWindow::render() {
  if (visible_) {
    if (ImGui::Begin("Pump Setup", &visible_)) {
      voltageSlider_->render();
      freqSlider_->render();
      valveOnOff_->render();
      ImGui::End();
    }
  }
}

void PumpWindow::setPumps() {
  for (int i = 0; i < NUM_PUMPS; i++)
    pp_->setVoltage(i + 1, (int16_t)voltageSlider_->get(i));
}

void PumpWindow::setFreq() {
  pp_->setFreq(freqSlider_->get(0));
}

// (ON = TRUE = Close = LOW, OFF = FALSE = Open = high)
void PumpWindow::setValves() {
  for (int i = 0; i < NUM_PUMPS; i++) {
    if (pp_->valveState[i]) { // close/ON -> open/OFF
      valveOnOff_->setLabel(i, "Close\nValve " + std::to_string(i + 1));
      pp_->setValve(i + 1, false);
    } else { // open/OFF -> close/ON
      valveOnOff_->setLabel(i, "Open\nValve " + std::to_string(i + 1));
      pp_->setValve(i + 1, true);
    }
  }
}

} // namespace gui
