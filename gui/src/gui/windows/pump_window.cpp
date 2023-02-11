#include "gui/windows/pump_window.hpp"

namespace gui {

PumpWindow::PumpWindow(std::shared_ptr<Pump> pp) : pp_(pp) {
  voltageSlider_ = std::make_unique<Slider<int>>("Voltage", 0, 250, "%d V", NUM_PUMPS,
                                                 ImVec2(40, 200), false, [this]() { setPumps(); });
  freqSlider_ = std::make_unique<Slider<int>>("Frequency", 0, 800, "%d Hz Freq", 1, ImVec2(40, 200),
                                              true, [this]() { setFreq(); });
  valveOnOff_ = std::make_unique<Button>(
      std::vector<std::string>{"Open\nValve\n1", "Open\nValve\n2", "Open\nValve\n3",
                               "Open\nValve\n4"},
      [this]() { setValves(); }, NUM_PUMPS, ImVec2(40, 0));
  controlOnOff_ = std::make_unique<Button>("Control On", [this]() {
    if (controlOnOff_->get(0))
      controlFlag_ = !controlFlag_;
    controlFlag_ ? controlOnOff_->setLabel(0, "Control Off")
                 : controlOnOff_->setLabel(0, "Control On");
    resetBtn_ = std::make_unique<Button>("Reset Pump", [this]() { resetPump(); });
  });
}

PumpWindow::~PumpWindow() {
  voltageSlider_.reset();
  freqSlider_.reset();
  valveOnOff_.reset();
  controlOnOff_.reset();
  resetBtn_.reset();
}

void PumpWindow::render() {
  if (ImGui::IsKeyPressed(ImGuiKey_P))
    visible_ = !visible_;
  if (visible_) {
    if (ImGui::Begin("Pump Setup", &visible_)) {
      voltageSlider_->render();
      freqSlider_->render();
      valveOnOff_->render();
      controlOnOff_->render();
      resetBtn_->render();
      ImGui::End();
    }
  }
}

void PumpWindow::setPumps() {
  if (!controlFlag_)
    return;
  for (int i = 0; i < NUM_PUMPS; i++)
    pp_->setVoltage(i + 1, (int16_t)voltageSlider_->get(i));
}

void PumpWindow::setFreq() {
  if (!controlFlag_)
    return;
  pp_->setFreq(freqSlider_->get(0));
}

// (ON = TRUE = Close = LOW, OFF = FALSE = Open = high)
void PumpWindow::setValves() {
  if (!controlFlag_)
    return;
  for (int i = 0; i < NUM_PUMPS; i++) {
    if (valveOnOff_->get(i)) {
      if (pp_->valveState[i]) { // close/ON -> open/OFF
        valveOnOff_->setLabel(i, "Close\nValve\n" + std::to_string(i + 1));
        pp_->setValve(i + 1, false);
      } else { // open/OFF -> close/ON
        valveOnOff_->setLabel(i, "Open\nValve\n" + std::to_string(i + 1));
        pp_->setValve(i + 1, true);
      }
    }
  }
}

void PumpWindow::resetPump() {
  if (resetBtn_->get(0)) {
    for (int i = 0; i < NUM_PUMPS; i++) {
      pp_->setValve(i + 1, false);
      valveOnOff_->setLabel(i, "Open\nValve\n" + std::to_string(i + 1));
      pp_->setVoltage(i + 1, 0);
      voltageSlider_->set(i, 0);
    }
    pp_->setFreq(0);
    freqSlider_->set(0, 0);
  }
}

} // namespace gui
