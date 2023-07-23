#include "gui/windows/pump_window.hpp"

namespace gui {

PumpWindow::PumpWindow(Pump *pp) : pp_(pp) {
  pumpUnits_ = pp_->getPumpType() == "FLUIGENT" ? "%d mbar" : "%d V";
  maxOutputSlider_ =
      std::make_unique<Slider<int>>("Max\nOutput", 30, 1000, nullptr, pumpUnits_, ImVec2(0, 0),
                                    true, [this]() { setMaxOutput(); });
  outputSlider_ = std::make_unique<SliderArray<float>>(
      "Pump\nOutput", 0, 250, &pp_->outputs, pp_->getNumPumps(), pumpUnits_, ImVec2(40, 200), false,
      [this]() { setPumps(); });
  syncToggle_ = std::make_unique<Toggle>("Sync pump 1 & 2", nullptr, [this]() { syncPumps(); });
  freqSlider_ = std::make_unique<Slider<int>>("Frequency", 0, 800, &pp_->freq, "%d Hz",
                                              ImVec2(0, 0), true, [this]() { setFreq(); });
  valveToggle_ = std::make_unique<ToggleArray>(
      "Valve On/Off", pp_->valveState, pp_->getNumPumps(), [this]() { setValves(); }, ImVec2(0, 0));
  controlToggle_ = std::make_unique<Toggle>("Control On/Off");
  resetBtn_ = std::make_unique<Button>("Reset Pump", nullptr, [this]() { resetPump(); });
}

PumpWindow::~PumpWindow() {
  maxOutputSlider_.reset();
  outputSlider_.reset();
  syncToggle_.reset();
  freqSlider_.reset();
  valveToggle_.reset();
  controlToggle_.reset();
  resetBtn_.reset();
}

void PumpWindow::render() {
  if (ImGui::GetIO().KeyAlt && ImGui::IsKeyPressed(ImGuiKey_P))
    visible_ = !visible_;
  if (visible_ && ImGui::Begin("Pump Setup", &visible_)) {
    maxOutputSlider_->render();
    outputSlider_->render();
    if (controlToggle_->get())
      for (int p = 0; p < pp_->getNumPumps(); ++p)
        ImGui::InputFloat(("P" + std::to_string(p + 1)).c_str(), &pp_->outputs[p], 0.1, 1);
    syncToggle_->render();
    if (pp_->getPumpType() == "BARTELS") {
      freqSlider_->render();
      valveToggle_->render();
    }
    controlToggle_->render();
    resetBtn_->render();
    ImGui::End();
  }
}

void PumpWindow::setPumps() {
  if (!controlToggle_->get())
    return;
  for (int i = 0; i < pp_->getNumPumps(); ++i)
    pp_->setOutput(i, outputSlider_->get(i));
}

void PumpWindow::syncPumps() { pp_->outputs[1] = pp_->outputs[0]; }

void PumpWindow::setFreq() {
  if (!controlToggle_->get())
    return;
  pp_->setFreq(freqSlider_->get());
}

// (ON = TRUE = Close = LOW, OFF = FALSE = Open = high)
void PumpWindow::setValves() {
  if (!controlToggle_->get())
    return;
  for (int i = 0; i < pp_->getNumPumps(); ++i)
    if (valveToggle_->changed(i))
      pp_->setValve(i, valveToggle_->get(i));
}

void PumpWindow::setMaxOutput() { outputSlider_->setMax(maxOutputSlider_->get()); }

void PumpWindow::resetPump() {
  if (resetBtn_->get()) {
    for (int i = 0; i < pp_->getNumPumps(); ++i) {
      if (pp_->getPumpType() == "BARTELS")
        pp_->setValve(i, false);
      pp_->setOutput(i, 0);
    }
    if (pp_->getPumpType() == "BARTELS")
      pp_->setFreq(0);
  }
}

} // namespace gui
