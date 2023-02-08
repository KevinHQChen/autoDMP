#pragma once

#include "ctrl/supervisor.hpp"
#include "gui/components/button.hpp"
#include "gui/components/checkboxarray.hpp"
#include "gui/components/dropdown.hpp"
#include "gui/components/sliderfloatarray.hpp"
#include "gui/components/sliderint.hpp"
#include "window.hpp"
#include <ctrl/state/state.hpp>

#include "gui/components/implot_helpers.hpp"
#include "implot.h"

#define NUM_CHANS 3

namespace gui {

class SysIdWindow : public Window {
  void previewExcitationSignal();
  void toggleExcitationSignal();
  void sendExcitationSignal();
  void stopExcitationSignal();
  void clearCtrlDataQueue();

  std::shared_ptr<Supervisor> sv_;

  std::unique_ptr<Button> excitationSignalPreviewBtn_;
  std::unique_ptr<Button> toggleExcitationSignalBtn_;
  std::unique_ptr<Button> sendExcitationSignalBtn_;
  std::unique_ptr<Button> stopExcitationSignalBtn_;
  std::unique_ptr<Button> clearDataBtn_;
  std::unique_ptr<Dropdown> excitationSignalDropdown_;
  std::unique_ptr<SliderFloatArray> minValSlider_;
  std::unique_ptr<SliderFloatArray> maxValSlider_;
  std::unique_ptr<SliderFloatArray> urefSlider_;
  std::unique_ptr<CheckboxArray> chSelect_;
  std::unique_ptr<SliderInt> numSampleSlider_;

  std::vector<std::string> excitationSignalTypes_ = {"sine", "square", "triangle", "sawtooth"};

  bool sysIDWindowVisible_ = false;
  float guiTime{0.0f}, history{30.0f};
  ScrollingBuffer u0, u1, u2, du0, du1, du2;
  ScrollingBuffer y0, y1, y2, yref0, yref1, yref2;
  std::vector<std::pair<ScrollingBuffer *, std::string>> sysidCtrlVecs{
      std::make_pair(&u0, "u0"), std::make_pair(&u1, "u1"), std::make_pair(&u2, "u2")};
  std::vector<std::pair<ScrollingBuffer *, std::string>> sysidMeasVecs{
      std::make_pair(&y0, "y0"), std::make_pair(&y1, "y1"), std::make_pair(&y2, "y2")};

public:
  SysIdWindow(std::shared_ptr<Supervisor> sv);
  ~SysIdWindow();
  void render() override;

  void plotVector3d(const char *plotName, const char *xAx, const char *yAx, double yMin,
                    double yMax, std::vector<std::pair<ScrollingBuffer *, std::string>> &vecs) {
    if (ImPlot::BeginPlot(plotName, ImVec2(-1, 300))) {
      ImPlot::SetupAxes(xAx, yAx); //, implotFlags, implotFlags);
      ImPlot::SetupAxisLimits(ImAxis_X1, guiTime - history, guiTime, ImGuiCond_Always);
      ImPlot::SetupAxisLimits(ImAxis_Y1, yMin, yMax);
      for (auto &vec : vecs)
        ImPlot::PlotLine(vec.second.c_str(), &vec.first->Data[0].x, &vec.first->Data[0].y,
                         vec.first->Data.size(), vec.first->Offset, 2 * sizeof(float));
      ImPlot::EndPlot();
    }
  }
};

} // namespace gui
