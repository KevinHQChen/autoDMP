#include "gui/windows/sysid_window.hpp"

namespace gui {

SysIdWindow::SysIdWindow(std::shared_ptr<Supervisor> sv) : sv_(sv) {
  chSelect_ = std::make_unique<CheckboxArray>("Channel", NUM_CHANS);
  numSampleSlider_ = std::make_unique<Slider<int>>("Num Samples", 0, 4000, &numSamples_);
  excitationSignalDropdown_ = std::make_unique<Dropdown>(
      "Excitation Signal Type:", excitationSignalTypes_, [this]() { generateExcitationSignal(); });
  minValSlider_ =
      std::make_unique<SliderArray<float>>("Excitation Signal Min Value", -10.0f, 0.0f, &minVal_);
  maxValSlider_ =
      std::make_unique<SliderArray<float>>("Excitation Signal Max Value", 0.0f, 10.0f, &maxVal_);
  urefSlider_ =
      std::make_unique<SliderArray<float>>("Control Signal Setpoint (uref)", 0.0f, 200.0f, &uref_);
  sendExcitationSignalBtn_ = std::make_unique<Button>("Send Excitation Signal", nullptr,
                                                      [this]() { sendExcitationSignal(); });
  stopExcitationSignalBtn_ = std::make_unique<Button>("Stop Excitation Signal", nullptr,
                                                      [this]() { stopExcitationSignal(); });
  clearDataBtn_ =
      std::make_unique<Button>("Clear ctrlDataQueue", nullptr, [this]() { clearCtrlDataQueue(); });
  // This is how you can add callbacks to the window
  // registerCallback([this]() { excitationSignalDropdown_->render(); });
}

SysIdWindow::~SysIdWindow() {
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
      ImGui::Separator();

      chSelect_->render();
      minValSlider_->render();
      maxValSlider_->render();
      urefSlider_->render();
      numSampleSlider_->render();
      excitationSignalDropdown_->render();
      sendExcitationSignalBtn_->render();
      stopExcitationSignalBtn_->render();
      clearDataBtn_->render();

      ImGui::End();
    }
  }

  if (sysIDWindowVisible_) {
    std::vector<float> uref = urefSlider_->get();

    sv_->startSysIDThread(Eigen::Vector3d(uref[0], uref[1], uref[2]), chSelect_->get(),
                          minValSlider_->get(), maxValSlider_->get(), excitationSignal_);

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
    }
  } else
    sv_->stopSysIDThread();
}

void SysIdWindow::generateExcitationSignal() {
  if (excitationSignalDropdown_->getItem() == "prbs") {
    info("Generating excitation signal...");
    py::gil_scoped_acquire acquire;
    py::eval_file("ctrl/scripts/sysid.py"); // import sysid functions
    py::object prbs = py::module::import("prbs").attr("prbs");
    py::list chSel;
    for (int i = 0; i < NUM_CHANS; i++)
      if (chSelect_->get()[i])
        chSel.append(chSelect_->get()[i]);

    excitationSignal_ =
        prbs(chSel, minValSlider_->get(), maxValSlider_->get(), numSampleSlider_->get())
            .cast<Eigen::MatrixXd>();
  }

  info("Excitation signal dimensions: {}x{}", excitationSignal_.rows(), excitationSignal_.cols());
}

void SysIdWindow::sendExcitationSignal() { sysIDWindowVisible_ = true; }

void SysIdWindow::stopExcitationSignal() { sysIDWindowVisible_ = false; }

void SysIdWindow::clearCtrlDataQueue() { sv_->ctrlDataQueuePtr->clearFile(); }

} // namespace gui
