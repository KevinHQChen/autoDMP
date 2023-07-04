#pragma once

#include "gui/components/button.hpp"
#include "gui/components/implot_helpers.hpp"
#include "gui/components/slider.hpp"
#include "window.hpp"

#include "ctrl/supervisor.hpp"

namespace gui {

class PlotWindow : public Window {
  Supervisor *sv_;
  ImProc *imProc_;

  float guiTime{0.0f}, history{30.0f};
  bool pausePlot{false};
  int plotHeight{200};

  int openAction = -1;

  std::vector<double> y_;
  std::vector<ScrollingBuffer *> u, uref, y, yhat, yref, ywt, theta, currTraj;

  void initDataVecs();
  void destroyDataVecs();

public:
  PlotWindow(ImProc *imProc, Supervisor *sv);
  ~PlotWindow();
  void render() override;
};

} // namespace gui
