#include "gui/windows/plot_window.hpp"

namespace gui {

PlotWindow::PlotWindow(ImProc *imProc, Supervisor *sv) : imProc_(imProc), sv_(sv) {}

PlotWindow::~PlotWindow() { destroyDataVecs(); }

void PlotWindow::render() {
  initDataVecs();
  if (ImGui::GetIO().KeyAlt && ImGui::IsKeyPressed(ImGuiKey_V))
    visible_ = !visible_;
  if (visible_ && ImGui::Begin("Real Time Plot", &visible_)) {
    if (!pausePlot) {
      // update all plotted data
      guiTime += ImGui::GetIO().DeltaTime;
      if (imProc_->started()) {
        for (int i = 0; i < 2 * imProc_->impConf.getNumChs(); ++i) {
          y_ = imProc_->getY();
          y[i]->AddPoint(guiTime, y_[i]);
        }
      }
      if (sv_->started()) {
        for (int i = 0; i < sv_->no; ++i)
          u[i]->AddPoint(guiTime, sv_->supOut.u[i]);
        for (int i = 0; i < 2 * sv_->no; ++i) {
          yhat[i]->AddPoint(guiTime, sv_->supOut.yhat[i]);
          ywt[i]->AddPoint(guiTime, sv_->supOut.ywt[i]);
          currTraj[i]->AddPoint(guiTime, sv_->supOut.currTraj[i]);
        }
        int np = 4;
        for (int i = 0; i < np * 2 * sv_->no; ++i)
          theta[i]->AddPoint(guiTime, sv_->supOut.theta[i]);
      }
    }

    openAction = -1;
    if (openAction != 1 || openAction == -1) {
      if (ImGui::Button("Open all"))
        openAction = 1;
    } else if (openAction != 0) {
      if (ImGui::Button("Close all"))
        openAction = 0;
    }
    // disable tree node indentation
    ImGui::PushStyleVar(ImGuiStyleVar_IndentSpacing, 0.0f);
    ImGui::Separator();

    if (openAction != -1)
      ImGui::SetNextItemOpen(openAction != 0);
    if (ImGui::TreeNode("Settings")) {
      ImGui::SliderFloat("History", &history, 1, 60, "%.1f s");
      ImGui::SliderInt("Plot Height", &plotHeight, 10, 300);
      if (!pausePlot) {
        if (ImGui::Button("Pause Plot"))
          pausePlot = true;
      } else {
        if (ImGui::Button("Resume Plot"))
          pausePlot = false;
      }
      ImGui::SameLine();
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
        for (auto &vec : currTraj)
          vec->Erase();
        for (auto &vec : theta)
          vec->Erase();
      }
      ImGui::TreePop();
    }
    ImGui::Separator();

    if (openAction != -1)
      ImGui::SetNextItemOpen(openAction != 0);
    if (ImGui::TreeNode("Output")) {
      plotVectorNN("##Output", "time (s)", "position (px)", -100, 100, {"y", "yhat", "currTraj"},
                   {y, yhat, currTraj}, plotHeight, guiTime, history);
      ImGui::TreePop();
    }
    ImGui::Separator();

    if (openAction != -1)
      ImGui::SetNextItemOpen(openAction != 0);
    if (ImGui::TreeNode("Input")) {
      plotVectorN("##Input", "time (s)", "input (mbar)", 0, 80, "u", u, plotHeight, guiTime,
                  history);
      ImGui::TreePop();
    }
    ImGui::Separator();

    if (openAction != -1)
      ImGui::SetNextItemOpen(openAction != 0);
    if (ImGui::TreeNode("Output Weight")) {
      plotVectorN("##Output Weight", "time (s)", "output weight", -0.1, 1.1, "ywt", ywt, plotHeight,
                  guiTime, history);
      ImGui::TreePop();
    }
    ImGui::Separator();

    if (openAction != -1)
      ImGui::SetNextItemOpen(openAction != 0);
    if (ImGui::TreeNode("Param Estimates")) {
      plotVectorN("##Param Estimates", "time (s)", "value", -1, 1, "th", theta, plotHeight, guiTime,
                  history);
      ImGui::TreePop();
    }
    ImGui::Separator();

    ImGui::PopStyleVar();
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
    if (!currTraj.empty())
      currTraj.clear();
  }
  if (sv_->started()) {
    if (yhat.empty())
      for (int ch = 0; ch < 2 * sv_->no; ++ch)
        yhat.push_back(new ScrollingBuffer());
    if (u.empty())
      for (int ch = 0; ch < sv_->no; ++ch)
        u.push_back(new ScrollingBuffer());
    if (ywt.empty())
      for (int ch = 0; ch < 2 * sv_->no; ++ch)
        ywt.push_back(new ScrollingBuffer());
    if (currTraj.empty())
      for (int ch = 0; ch < 2 * sv_->no; ++ch)
        currTraj.push_back(new ScrollingBuffer());
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
  for (auto &vec : currTraj)
    delete vec;
  for (auto &vec : theta)
    delete vec;
}

} // namespace gui
