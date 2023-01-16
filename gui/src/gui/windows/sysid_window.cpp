#include "gui/windows/sysid_window.hpp"

namespace gui {

SysIdWindow::SysIdWindow(std::shared_ptr<Supervisor> sv) : sv_(sv) {}

SysIdWindow::~SysIdWindow() {}

void SysIdWindow::render() {
  if (visible_) {
    sv_->startSysIDThread();

    if (ImGui::Begin("SysID", &visible_)) {
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

} // namespace gui
