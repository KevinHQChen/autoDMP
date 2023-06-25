#include "gui/windows/plot_window.hpp"

namespace gui {

PlotWindow::PlotWindow(ImProc *imProc, Supervisor *sv) : imProc_(imProc), sv_(sv) {}

PlotWindow::~PlotWindow() { destroyDataVecs(); }

void PlotWindow::render() {
  initDataVecs();
  if (ImGui::IsKeyPressed(ImGuiKey_V))
    visible_ = !visible_;
  if (visible_ && ImGui::Begin("Real Time Plot", &visible_)) {
    if (!pausePlot) {
      guiTime += ImGui::GetIO().DeltaTime;
      if (imProc_->started()) {
        for (int i = 0; i < 2 * imProc_->impConf.getNumChs(); ++i) {
          y_ = imProc_->getY();
          y[i]->AddPoint(guiTime, y_[i]);
        }
      }

      if (sv_->started()) {
        for (int i = 0; i < 2 * sv_->no; ++i) {
          yhat[i]->AddPoint(guiTime, sv_->supOut.yhat[i]);
          u[i]->AddPoint(guiTime, sv_->supOut.u[i]);
          ywt[i]->AddPoint(guiTime, sv_->supOut.ywt[i]);
        }
        int np = 4;
        for (int i = 0; i < np * 2 * sv_->no; ++i)
          theta[i]->AddPoint(guiTime, sv_->supOut.theta[i]);
      }
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

    plotVectorNN("Output", "time (s)", "position (px)", 0, 80, {"y", "yhat"}, {y, yhat}, guiTime,
                 history);
    plotVectorN("Input", "time (s)", "input (mbar)", 0, 80, "u", u, guiTime, history);
    plotVectorN("Output Weight", "time (s)", "output weight", 0, 1, "ywt", ywt, guiTime, history);
    plotVectorN("Param Estimates", "time (s)", "value", -1, 1, "th", theta, guiTime, history);
    ImGui::End();
  }
}

void PlotWindow::initDataVecs() {
  if (!imProc_->started() && !y.empty())
    y.clear();
  if (imProc_->started() && y.empty())
    for (int ch = 0; ch < 2 * imProc_->impConf.getNumChs(); ++ch)
      y.push_back(new ScrollingBuffer());

  if (!sv_->started()) {
    if (!yhat.empty())
      yhat.clear();
    if (!u.empty())
      u.clear();
    if (!ywt.empty())
      ywt.clear();
    if (!theta.empty())
      theta.clear();
  }
  if (sv_->started()) {
    if (yhat.empty())
      for (int ch = 0; ch < 2 * sv_->no; ++ch)
        yhat.push_back(new ScrollingBuffer());
    if (u.empty())
      for (int ch = 0; ch < 2 * sv_->no; ++ch)
        u.push_back(new ScrollingBuffer());
    if (ywt.empty())
      for (int ch = 0; ch < 2 * sv_->no; ++ch)
        ywt.push_back(new ScrollingBuffer());
    int np = 4;
    if (theta.empty())
      for (int i = 0; i < np * 2 * sv_->no; ++i)
        theta.push_back(new ScrollingBuffer());
  }
}

void PlotWindow::destroyDataVecs() {
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

} // namespace gui
