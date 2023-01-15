#include "gui/windows/sysid_setup_window.hpp"

namespace gui {

SysIdSetupWindow::SysIdSetupWindow() {
  excitationSignalPreviewButton_ = std::make_unique<Button>(
      "Preview Excitation Signal", [this]() { previewExcitationSignal(); });
  registerCallback([this]() { excitationSignalPreviewButton_->render(); });

  toggleExcitationSignalButton_ =
      std::make_unique<Button>("Toggle Excitation Signal", [this]() { toggleExcitationSignal(); });
  registerCallback([this]() { toggleExcitationSignalButton_->render(); });

  sendExcitationSignalButton_ =
      std::make_unique<Button>("Send Excitation Signal", [this]() { sendExcitationSignal(); });
  registerCallback([this]() { sendExcitationSignalButton_->render(); });
}

SysIdSetupWindow::~SysIdSetupWindow() {
  excitationSignalPreviewButton_.reset();
  toggleExcitationSignalButton_.reset();
  sendExcitationSignalButton_.reset();
  excitationSignalDropdown_.reset();
}

void SysIdSetupWindow::render() {
  ImGui::Begin("SysID Setup");
  ImGui::Text("Configure Excitation Signal");

  excitationSignalDropdown_->render();
  excitationSignalPreviewButton_->render();
  toggleExcitationSignalButton_->render();

  ImGui::Text("Excitation Signal Max/Min Values:");
  ImGui::BeginGroup();
  for (int i = 0; i < 3; i++) {
    std::string maxLabel = "Channel " + std::to_string(i) + " Max";
    std::string minLabel = "Channel " + std::to_string(i) + " Min";
    ImGui::SliderFloat(maxLabel.c_str(), &excitationSignalMaxValues_[i], 0.0f, 1.0f);
    ImGui::SliderFloat(minLabel.c_str(), &excitationSignalMinValues_[i], 0.0f, 1.0f);
  }
  ImGui::EndGroup();

  ImGui::Text("Control Signal Setpoint (uref):");
  ImGui::BeginGroup();
  for (int i = 0; i < 3; i++) {
    std::string label = "Channel " + std::to_string(i);
    ImGui::SliderFloat(label.c_str(), &controlSignalSetpoints_[i], 0.0f, 1.0f);
  }
  ImGui::EndGroup();

  sendExcitationSignalButton_->render();

  ImGui::End();
}

void SysIdSetupWindow::previewExcitationSignal() { info("Previewing excitation signal"); }

void SysIdSetupWindow::toggleExcitationSignal() { info("Toggling excitation signal"); }

void SysIdSetupWindow::sendExcitationSignal() { info("Sending excitation signal"); }

} // namespace gui
