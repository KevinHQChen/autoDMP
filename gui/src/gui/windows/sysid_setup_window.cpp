#include "gui/windows/sysid_setup_window.hpp"

namespace gui {

SysIdSetupWindow::SysIdSetupWindow() {
  excitationSignalDropdown_ =
      std::make_unique<Dropdown>("Excitation Signal Type:", excitationSignalTypes_,
                                 [this]() { info("Excitation Signal Type"); });
  excitationSignalPreviewButton_ = std::make_unique<Button>(
      "Preview Excitation Signal", [this]() { previewExcitationSignal(); });
  toggleExcitationSignalButton_ =
      std::make_unique<Button>("Toggle Excitation Signal", [this]() { toggleExcitationSignal(); });
  minValSlider_ = std::make_unique<SliderFloatArray>(
      "Excitation Signal Min Value", 3, 0.0f, 1.0f, [this]() { info("Min Value"); });
  maxValSlider_ = std::make_unique<SliderFloatArray>(
      "Excitation Signal Max Value", 3, 0.0f, 1.0f, [this]() { info("Max Value"); });
  urefSlider_ = std::make_unique<SliderFloatArray>(
      "Control Signal Setpoint (uref)", 3, 0.0f, 1.0f, [this]() { info("Uref"); });
  sendExcitationSignalButton_ =
      std::make_unique<Button>("Send Excitation Signal", [this]() { sendExcitationSignal(); });
  // This is how you can add callbacks to the window
  // registerCallback([this]() { excitationSignalDropdown_->render(); });
}

SysIdSetupWindow::~SysIdSetupWindow() {
  excitationSignalPreviewButton_.reset();
  toggleExcitationSignalButton_.reset();
  sendExcitationSignalButton_.reset();
  excitationSignalDropdown_.reset();
  minValSlider_.reset();
  maxValSlider_.reset();
  urefSlider_.reset();
}

void SysIdSetupWindow::render() {
  ImGui::Begin("SysID Setup");
  ImGui::Text("Configure Excitation Signal");
  excitationSignalDropdown_->render();
  excitationSignalPreviewButton_->render();
  toggleExcitationSignalButton_->render();
  minValSlider_->render();
  maxValSlider_->render();
  urefSlider_->render();
  sendExcitationSignalButton_->render();
  ImGui::End();
}

void SysIdSetupWindow::previewExcitationSignal() { info("Previewing excitation signal"); }

void SysIdSetupWindow::toggleExcitationSignal() { info("Toggling excitation signal"); }

void SysIdSetupWindow::sendExcitationSignal() { info("Sending excitation signal"); }

} // namespace gui
