#pragma once

#include "ctrl/supervisor.hpp"
#include "gui/components/button.hpp"
#include "gui/components/checkboxarray.hpp"
#include "gui/components/dropdown.hpp"
#include "gui/components/slider.hpp"
#include "window.hpp"
#include <ctrl/state/state.hpp>

#include "gui/components/implot_helpers.hpp"
#include "implot/implot.h"
#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace gui {

namespace py = pybind11;
using namespace py::literals;

class SysIdWindow : public Window {
  void generateExcitationSignal();
  void sendExcitationSignal();
  void stopExcitationSignal();
  void clearCtrlDataQueue();

  Supervisor *sv_;

  std::unique_ptr<Button> sendExcitationSignalBtn_;
  std::unique_ptr<Button> stopExcitationSignalBtn_;
  std::unique_ptr<Button> clearDataBtn_;
  std::unique_ptr<Dropdown> excitationSignalDropdown_;
  std::unique_ptr<SliderArray<float>> minValSlider_;
  std::unique_ptr<SliderArray<float>> maxValSlider_;
  std::unique_ptr<SliderArray<float>> urefSlider_;
  std::unique_ptr<CheckboxArray> chSelect_;
  std::unique_ptr<Slider<int>> numSampleSlider_;
  std::unique_ptr<Toggle> flipSignalToggle_;

  std::vector<std::string> excitationSignalTypes_ = {"prbs", "sine", "square", "triangle",
                                                     "sawtooth"};
  Eigen::MatrixXd excitationSignal_;
  std::vector<double> timeVec_, u0Vec_, u1Vec_, u2Vec_;
  // py::object prbs;

  int numSamples_ = 1000;
  std::vector<float> minVal_{std::vector<float>(NUM_CHANS, 0.0f)};
  std::vector<float> maxVal_{std::vector<float>(NUM_CHANS, 1.0f)};
  std::vector<float> uref_{std::vector<float>(NUM_CHANS, 1.0f)};

  float guiTime{0.0f}, history{30.0f};
  ScrollingBuffer u0, u1, u2, du0, du1, du2;
  ScrollingBuffer y0, y1, y2, yref0, yref1, yref2;
  std::vector<std::pair<ScrollingBuffer *, std::string>> sysidCtrlVecs{
      std::make_pair(&u0, "u0"), std::make_pair(&u1, "u1"), std::make_pair(&u2, "u2")};
  std::vector<std::pair<ScrollingBuffer *, std::string>> sysidMeasVecs{
      std::make_pair(&y0, "y0"), std::make_pair(&y1, "y1"), std::make_pair(&y2, "y2")};

public:
  bool sysIDWindowVisible_{false};

  SysIdWindow(Supervisor *sv);
  ~SysIdWindow();
  void render() override;
};

} // namespace gui
