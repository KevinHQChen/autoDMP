#include "gui/windows/pump_window.hpp"

namespace gui {

PumpWindow::PumpWindow(Pump *pp) : pp_(pp) {
  maxVoltageSlider_ = std::make_unique<Slider<int>>(
      "Max\nVoltage", 30, 1000, nullptr, "%d V", ImVec2(0, 0), true, [this]() { setMaxVoltage(); });
  voltageSlider_ =
      std::make_unique<SliderArray<int>>("Pump\nVoltage", 0, 250, &pp_->pumpVoltages, NUM_PUMPS,
                                         "%d V", ImVec2(40, 200), false, [this]() { setPumps(); });
  syncToggle_ = std::make_unique<Toggle>("Sync pump 1 & 2", nullptr, [this]() { syncPumps(); });
  freqSlider_ = std::make_unique<Slider<int>>("Frequency", 0, 800, &pp_->freq, "%d Hz",
                                              ImVec2(0, 0), true, [this]() { setFreq(); });
  valveToggle_ = std::make_unique<ToggleArray>(
      "Valve On/Off", pp_->valveState, NUM_PUMPS, [this]() { setValves(); }, ImVec2(0, 0));
  controlToggle_ = std::make_unique<Toggle>("Control On/Off");
  resetBtn_ = std::make_unique<Button>("Reset Pump", nullptr, [this]() { resetPump(); });
}

PumpWindow::~PumpWindow() {
  maxVoltageSlider_.reset();
  voltageSlider_.reset();
  syncToggle_.reset();
  freqSlider_.reset();
  valveToggle_.reset();
  controlToggle_.reset();
  resetBtn_.reset();
}

void PumpWindow::render() {
  if (ImGui::IsKeyPressed(ImGuiKey_P))
    visible_ = !visible_;
  if (visible_) {
    if (ImGui::Begin("Pump Setup", &visible_)) {
      maxVoltageSlider_->render();
      voltageSlider_->render();
      syncToggle_->render();
      freqSlider_->render();
      valveToggle_->render();
      controlToggle_->render();
      resetBtn_->render();
      ImGui::End();
    }
  }
}

void PumpWindow::setPumps() {
  if (!controlToggle_->get())
    return;
  for (int i = 0; i < NUM_PUMPS; i++)
    pp_->setVoltage(i, (int16_t)voltageSlider_->get(i));
}

void PumpWindow::syncPumps() {
  if (syncToggle_->get())
    pp_->pumpVoltages[1] = pp_->pumpVoltages[0];
}

void PumpWindow::setFreq() {
  if (!controlToggle_->get())
    return;
  pp_->setFreq(freqSlider_->get());
}

// (ON = TRUE = Close = LOW, OFF = FALSE = Open = high)
void PumpWindow::setValves() {
  if (!controlToggle_->get())
    return;
  for (int i = 0; i < NUM_PUMPS; i++)
    if (valveToggle_->changed(i))
      pp_->setValve(i, valveToggle_->get(i));
}

void PumpWindow::setMaxVoltage() { voltageSlider_->setMax(maxVoltageSlider_->get()); }

void PumpWindow::resetPump() {
  if (resetBtn_->get()) {
    for (int i = 0; i < NUM_PUMPS; i++) {
      pp_->setValve(i, false);
      pp_->setVoltage(i, 0);
    }
    pp_->setFreq(0);
  }
}

} // namespace gui
