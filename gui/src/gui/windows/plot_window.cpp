#include "gui/windows/plot_window.hpp"

namespace gui {

PlotWindow::PlotWindow(Supervisor *sv) : sv_(sv) {}

PlotWindow::~PlotWindow() {}

void PlotWindow::render() {
  if (ImGui::IsKeyPressed(ImGuiKey_V))
    visible_ = !visible_;
  if (visible_) {
    if (ImGui::Begin("Real Time Plot", &visible_)) {
      if (!pausePlot) {
        guiTime += ImGui::GetIO().DeltaTime;
        u0.AddPoint(guiTime, sv_->supOut.u[0]);
        u1.AddPoint(guiTime, sv_->supOut.u[1]);
        u2.AddPoint(guiTime, sv_->supOut.u[2]);
        uref0.AddPoint(guiTime, sv_->supOut.uref[0]);
        uref1.AddPoint(guiTime, sv_->supOut.uref[1]);
        uref2.AddPoint(guiTime, sv_->supOut.uref[2]);

        y0.AddPoint(guiTime, sv_->supIn.ymeas[0]);
        y1.AddPoint(guiTime, sv_->supIn.ymeas[1]);
        y2.AddPoint(guiTime, sv_->supIn.ymeas[2]);
        yhat0.AddPoint(guiTime, sv_->supOut.yhat[0]);
        yhat1.AddPoint(guiTime, sv_->supOut.yhat[1]);
        yhat2.AddPoint(guiTime, sv_->supOut.yhat[2]);
        yref0.AddPoint(guiTime, sv_->supOut.currTraj[0]);
        yref1.AddPoint(guiTime, sv_->supOut.currTraj[1]);
        yref2.AddPoint(guiTime, sv_->supOut.currTraj[2]);

        b11.AddPoint(guiTime, sv_->supOut.B_a[0]);
        b21.AddPoint(guiTime, sv_->supOut.B_a[1]);
        b31.AddPoint(guiTime, sv_->supOut.B_a[2]);
        b12.AddPoint(guiTime, sv_->supOut.B_a[3]);
        b22.AddPoint(guiTime, sv_->supOut.B_a[4]);
        b32.AddPoint(guiTime, sv_->supOut.B_a[5]);
        b13.AddPoint(guiTime, sv_->supOut.B_a[6]);
        b23.AddPoint(guiTime, sv_->supOut.B_a[7]);
        b33.AddPoint(guiTime, sv_->supOut.B_a[8]);
      }

      ImGui::SliderFloat("History", &history, 1, 60, "%.1f s");
      if (!pausePlot)
        if (ImGui::Button("Pause Plot"))
          pausePlot = true;
      if (pausePlot)
        if (ImGui::Button("Resume Plot"))
          pausePlot = false;
      if (ImGui::Button("Erase Data")) {
        guiTime = 0;
        for (auto &vec : std::vector<ScrollingBuffer>{u0,    u1,    u2,    y0,    y1,    y2,  yhat0,
                                                      yhat1, yhat2, yref0, yref1, yref2, b11, b21,
                                                      b31,   b12,   b22,   b32,   b13,   b23, b33})
          vec.Erase();
      }

      plotVector3d("##Control Input", "time (s)", "voltage (V)", 0, 90, ctrlVecs, guiTime,
                   history);
      plotVector3d("##Measured Output", "time (s)", "position (px)", 0, 750, measVecs, guiTime,
                   history);
      plotVector3d("##Param Estimates", "time (s)", "value", -1, 1, paramVecs, guiTime,
                   history);
      ImGui::End();
    }
  }
}

} // namespace gui
