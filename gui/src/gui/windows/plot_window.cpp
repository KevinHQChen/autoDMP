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

        y0.AddPoint(guiTime, sv_->supIn.y[0]);
        y1.AddPoint(guiTime, sv_->supIn.y[1]);
        y2.AddPoint(guiTime, sv_->supIn.y[2]);
        yhat0.AddPoint(guiTime, sv_->supOut.yhat[0]);
        yhat1.AddPoint(guiTime, sv_->supOut.yhat[1]);
        yhat2.AddPoint(guiTime, sv_->supOut.yhat[2]);
        yref0.AddPoint(guiTime, sv_->supOut.currTraj[0]);
        yref1.AddPoint(guiTime, sv_->supOut.currTraj[1]);
        yref2.AddPoint(guiTime, sv_->supOut.currTraj[2]);
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
