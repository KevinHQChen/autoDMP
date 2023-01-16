#include "gui/windows/sysid_window.hpp"

namespace gui {

SysIdWindow::SysIdWindow(std::shared_ptr<Supervisor> sv) : sv_(sv) {
  chSelect_ = std::make_unique<CheckboxArray>("Channel", NUM_CHANS);
  numSampleSlider_ = std::make_unique<SliderInt>("Num Samples", 0, 4000);
  excitationSignalDropdown_ =
      std::make_unique<Dropdown>("Excitation Signal Type:", excitationSignalTypes_);
  excitationSignalPreviewBtn_ = std::make_unique<Button>("Preview Excitation Signal",
                                                         [this]() { previewExcitationSignal(); });
  toggleExcitationSignalBtn_ =
      std::make_unique<Button>("Toggle Excitation Signal", [this]() { toggleExcitationSignal(); });
  minValSlider_ = std::make_unique<SliderFloatArray>("Excitation Signal Min Value", 3, 0.0f, 1.0f);
  maxValSlider_ = std::make_unique<SliderFloatArray>("Excitation Signal Max Value", 3, 0.0f, 1.0f);
  urefSlider_ = std::make_unique<SliderFloatArray>("Control Signal Setpoint (uref)", 3, 0.0f, 1.0f);
  sendExcitationSignalBtn_ =
      std::make_unique<Button>("Send Excitation Signal", [this]() { sendExcitationSignal(); });
  stopExcitationSignalBtn_ =
      std::make_unique<Button>("Stop Excitation Signal", [this]() { stopExcitationSignal(); });
  clearDataBtn_ =
      std::make_unique<Button>("Clear ctrlDataQueue", [this]() { clearCtrlDataQueue(); });
  // This is how you can add callbacks to the window
  // registerCallback([this]() { excitationSignalDropdown_->render(); });
}

SysIdWindow::~SysIdWindow() {
  excitationSignalPreviewBtn_.reset();
  toggleExcitationSignalBtn_.reset();
  sendExcitationSignalBtn_.reset();
  excitationSignalDropdown_.reset();
  minValSlider_.reset();
  maxValSlider_.reset();
  urefSlider_.reset();
  stopExcitationSignalBtn_.reset();
  clearDataBtn_.reset();
  chSelect_.reset();
  numSampleSlider_.reset();
}

void SysIdWindow::render() {
  if (visible_) {
    if (ImGui::Begin("SysId Setup", &visible_)) {
      ImGui::Text("Configure Excitation Signal");

      chSelect_->render();
      minValSlider_->render();
      maxValSlider_->render();
      urefSlider_->render();
      numSampleSlider_->render();
      excitationSignalDropdown_->render();
      excitationSignalPreviewBtn_->render();
      toggleExcitationSignalBtn_->render();
      sendExcitationSignalBtn_->render();
      stopExcitationSignalBtn_->render();
      clearDataBtn_->render();

      ImGui::End();
    }
  }

  if (sysIDWindowVisible_) {
    float *uref = urefSlider_->getValues();

    sv_->startSysIDThread(Eigen::Vector3d(uref[0], uref[1], uref[2]), chSelect_->get(),
                          minValSlider_->getValues(), maxValSlider_->getValues(),
                          numSampleSlider_->get());

    if (ImGui::Begin("SysID", &sysIDWindowVisible_)) {
      guiTime += ImGui::GetIO().DeltaTime;
      u0.AddPoint(guiTime, sv_->currState_->u(0));
      u1.AddPoint(guiTime, sv_->currState_->u(1));
      u2.AddPoint(guiTime, sv_->currState_->u(2));
      y0.AddPoint(guiTime, sv_->currState_->y(0));
      y1.AddPoint(guiTime, sv_->currState_->y(1));
      y2.AddPoint(guiTime, sv_->currState_->y(2));

      ImGui::SliderFloat("History", &history, 1, 30, "%.1f s");

      plotVector3d("##Control Input", "time (s)", "voltage (V)", 0, 250, sysidCtrlVecs);
      plotVector3d("##Measured Output", "time (s)", "position (px)", -500, 500, sysidMeasVecs);
      ImGui::End();
    } else
      sv_->stopSysIDThread();
  }
}

void SysIdWindow::previewExcitationSignal() { info("Previewing excitation signal"); }

void SysIdWindow::toggleExcitationSignal() { info("Toggling excitation signal"); }

void SysIdWindow::sendExcitationSignal() { sysIDWindowVisible_ = true; }

void SysIdWindow::stopExcitationSignal() { sysIDWindowVisible_ = false; }

void SysIdWindow::clearCtrlDataQueue() { sv_->ctrlDataQueuePtr->clearFile(); }

} // namespace gui
