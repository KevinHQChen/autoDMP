#pragma once

#include "gui/components/button.hpp"
#include "gui/components/implot_helpers.hpp"
#include "gui/components/slider.hpp"
#include "window.hpp"

#include "ctrl/supervisor.hpp"

namespace gui {

class PlotWindow : public Window {
  Supervisor *sv_;

  float guiTime{0.0f}, history{30.0f};

  std::vector<ScrollingBuffer *> u, uref, y, yhat, yref, ywt, theta;

public:
  bool pausePlot{false};

  PlotWindow(Supervisor *sv);
  ~PlotWindow();
  void render() override;
};

} // namespace gui
