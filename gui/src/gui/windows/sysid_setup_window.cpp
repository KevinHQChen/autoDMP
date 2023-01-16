#include "gui/windows/sysid_setup_window.hpp"

namespace gui {

SysIdSetupWindow::SysIdSetupWindow(std::shared_ptr<SysIdWindow> sysIDWindow,
                                   std::shared_ptr<Supervisor> sv)
    : sv_(sv), sysIdWindow_(std::weak_ptr(sysIDWindow)) {
  excitationSignalDropdown_ =
      std::make_unique<Dropdown>("Excitation Signal Type:", excitationSignalTypes_,
                                 [this]() { info("Excitation Signal Type"); });
  excitationSignalPreviewBtn_ = std::make_unique<Button>("Preview Excitation Signal",
                                                         [this]() { previewExcitationSignal(); });
  toggleExcitationSignalBtn_ =
      std::make_unique<Button>("Toggle Excitation Signal", [this]() { toggleExcitationSignal(); });
  minValSlider_ = std::make_unique<SliderFloatArray>("Excitation Signal Min Value", 3, 0.0f, 1.0f,
                                                     [this]() { info("Min Value"); });
  maxValSlider_ = std::make_unique<SliderFloatArray>("Excitation Signal Max Value", 3, 0.0f, 1.0f,
                                                     [this]() { info("Max Value"); });
  urefSlider_ = std::make_unique<SliderFloatArray>("Control Signal Setpoint (uref)", 3, 0.0f, 1.0f,
                                                   [this]() { info("Uref"); });
  sendExcitationSignalBtn_ =
      std::make_unique<Button>("Send Excitation Signal", [this]() { sendExcitationSignal(); });
  stopExcitationSignalBtn_ =
      std::make_unique<Button>("Stop Excitation Signal", [this]() { stopExcitationSignal(); });
  clearDataBtn_ =
      std::make_unique<Button>("Clear ctrlDataQueue", [this]() { clearCtrlDataQueue(); });
  // This is how you can add callbacks to the window
  // registerCallback([this]() { excitationSignalDropdown_->render(); });
}

SysIdSetupWindow::~SysIdSetupWindow() {
  excitationSignalPreviewBtn_.reset();
  toggleExcitationSignalBtn_.reset();
  sendExcitationSignalBtn_.reset();
  excitationSignalDropdown_.reset();
  minValSlider_.reset();
  maxValSlider_.reset();
  urefSlider_.reset();
  stopExcitationSignalBtn_.reset();
  clearDataBtn_.reset();
}

void SysIdSetupWindow::render() {
  if (visible_) {
    if (ImGui::Begin("SysId Setup", &visible_)) {
      ImGui::Text("Configure Excitation Signal");

      excitationSignalDropdown_->render();
      excitationSignalPreviewBtn_->render();
      toggleExcitationSignalBtn_->render();
      minValSlider_->render();
      maxValSlider_->render();
      urefSlider_->render();
      sendExcitationSignalBtn_->render();
      stopExcitationSignalBtn_->render();
      clearDataBtn_->render();

      ImGui::End();
    }
  }
}

void SysIdSetupWindow::previewExcitationSignal() { info("Previewing excitation signal"); }

void SysIdSetupWindow::toggleExcitationSignal() { info("Toggling excitation signal"); }

void SysIdSetupWindow::sendExcitationSignal() {
  if (auto sharedSysIdWindow_ = sysIdWindow_.lock())
    sharedSysIdWindow_->visible_ = true;
}

void SysIdSetupWindow::stopExcitationSignal() {
  if (auto sharedSysIdWindow_ = sysIdWindow_.lock())
    sharedSysIdWindow_->visible_ = false;
}

void SysIdSetupWindow::clearCtrlDataQueue() {
  if (auto sharedSysIdWindow_ = sysIdWindow_.lock())
    sv_->ctrlDataQueuePtr->clearFile();
}

} // namespace gui
