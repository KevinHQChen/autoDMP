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

  ScrollingBuffer u0, u1, u2, uref0, uref1, uref2;
  std::vector<std::pair<ScrollingBuffer *, std::string>> ctrlVecs{
      std::make_pair(&u0, "u0"),       std::make_pair(&u1, "u1"),
      std::make_pair(&u2, "u2"),       std::make_pair(&uref0, "uref0"),
      std::make_pair(&uref1, "uref1"), std::make_pair(&uref2, "uref2")};

  ScrollingBuffer y0, y1, y2, yhat0, yhat1, yhat2, yref0, yref1, yref2;
  std::vector<std::pair<ScrollingBuffer *, std::string>> measVecs{
      std::make_pair(&y0, "y0"),       std::make_pair(&y1, "y1"),
      std::make_pair(&y2, "y2"),       std::make_pair(&yhat0, "yhat0"),
      std::make_pair(&yhat1, "yhat1"), std::make_pair(&yhat2, "yhat2"),
      std::make_pair(&yref0, "yref0"), std::make_pair(&yref1, "yref1"),
      std::make_pair(&yref2, "yref2")};

  ScrollingBuffer b11, b21, b31, b12, b22, b32, b13, b23, b33;
  std::vector<std::pair<ScrollingBuffer *, std::string>> paramVecs{
      std::make_pair(&b11, "param0"), std::make_pair(&b21, "param1"),
      std::make_pair(&b31, "param2"), std::make_pair(&b12, "param3"),
      std::make_pair(&b22, "param4"), std::make_pair(&b32, "param5"),
      std::make_pair(&b13, "param6"), std::make_pair(&b23, "param7"),
      std::make_pair(&b33, "param8")};

public:
  bool pausePlot{false};

  PlotWindow(Supervisor *sv);
  ~PlotWindow();
  void render() override;
};

} // namespace gui
