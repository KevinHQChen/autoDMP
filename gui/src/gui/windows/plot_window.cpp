#include "gui/windows/plot_window.hpp"

namespace gui {

PlotWindow::PlotWindow(Supervisor *sv) : sv_(sv) {
  for (int ch = 0; ch < 2 * NUM_CHANS; ++ch) {
    y.push_back(new ScrollingBuffer());
    yhat.push_back(new ScrollingBuffer());

    u.push_back(new ScrollingBuffer());

    ywt.push_back(new ScrollingBuffer());
  }

  int np = 4;
  for (int i = 0; i < np * 2 * NUM_CHANS; ++i)
    theta.push_back(new ScrollingBuffer());
}

PlotWindow::~PlotWindow() {
  for (auto &vec : y)
    delete vec;
  for (auto &vec : yhat)
    delete vec;
  for (auto &vec : u)
    delete vec;
  for (auto &vec : ywt)
    delete vec;
  for (auto &vec : theta)
    delete vec;
}

void PlotWindow::render() {
  if (ImGui::IsKeyPressed(ImGuiKey_V))
    visible_ = !visible_;
  if (visible_) {
    if (ImGui::Begin("Real Time Plot", &visible_)) {
      if (!pausePlot) {
        guiTime += ImGui::GetIO().DeltaTime;
        for (int ch = 0; ch < 2 * NUM_CHANS; ++ch) {
          y[ch]->AddPoint(guiTime, sv_->supIn.y[ch]);
          yhat[ch]->AddPoint(guiTime, sv_->supOut.yhat[ch]);
          u[ch]->AddPoint(guiTime, sv_->supOut.u[ch]);
          ywt[ch]->AddPoint(guiTime, sv_->supOut.ywt[ch]);
        }
        int np = 4;
        for (int i = 0; i < np * 2 * NUM_CHANS; ++i)
          theta[i]->AddPoint(guiTime, sv_->supOut.theta[i]);
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
        for (auto &vec : y)
          vec->Erase();
        for (auto &vec : yhat)
          vec->Erase();
        for (auto &vec : u)
          vec->Erase();
        for (auto &vec : ywt)
          vec->Erase();
        for (auto &vec : theta)
          vec->Erase();
      }

      plotVectorNN("##Output", "time (s)", "position (px)", 0, 80, {"y", "yhat"}, {y, yhat},
                   guiTime, history);
      plotVectorN("##Input", "time (s)", "input (mbar)", 0, 80, "u", u, guiTime, history);
      plotVectorN("##Output Weight", "time (s)", "output weight", 0, 1, "ywt", ywt, guiTime,
                  history);
      plotVectorN("##Param Estimates", "time (s)", "value", -1, 1, "th", theta, guiTime, history);
      ImGui::End();
    }
  }
}

} // namespace gui
