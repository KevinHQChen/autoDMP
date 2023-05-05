#pragma once

#include "gui/windows.hpp"
#include "imgui_md_wrapper.h"
#include "immapp/immapp.h"
#include "implot/implot.h"

#include "ctrl/supervisor.hpp"
#include "imcap/imcap.hpp"
#include "improc/improc.hpp"
#include "pump/pump.hpp"
#include "util/util.hpp"

#include <cmath>
#include <cstdio>

namespace py = pybind11;
using namespace py::literals;

class GUI {
  HelloImGui::RunnerParams runnerParams;
  ImmApp::AddOnsParams addOnsParams;

  ordered_value conf;
  guiConfig guiConf;

  ImCap *imCap_;
  ImProc *imProc_;
  Pump *pump_;
  Supervisor *sv_;

  std::shared_ptr<gui::SysIdWindow> sysIDWindow_;
  std::shared_ptr<gui::PumpWindow> pumpWindow_;
  std::shared_ptr<gui::ImProcWindow> imProcWindow_;
  std::shared_ptr<gui::CtrlWindow> ctrlWindow_;

  // real time plotting
  // ImPlotAxisFlags implotFlags = ImPlotAxisFlags_NoTickLabels;
  // float guiTime{0.0f}, history{30.0f};
  // ScrollingBuffer u0, u1, u2;
  // ScrollingBuffer y0, y1, y2, yhat0, yhat1, yhat2, yref0, yref1, yref2;
  // std::vector<std::pair<ScrollingBuffer *, std::string>> sysidCtrlVecs{
  //     std::make_pair(&u0, "u0"), std::make_pair(&u1, "u1"), std::make_pair(&u2, "u2")};
  // std::vector<std::pair<ScrollingBuffer *, std::string>> sysidMeasVecs{
  //     std::make_pair(&y0, "y0"), std::make_pair(&y1, "y1"), std::make_pair(&y2, "y2")};
  // std::vector<std::pair<ScrollingBuffer *, std::string>> ctrlVecs{
  //     std::make_pair(&u0, "u0"), std::make_pair(&u1, "u1"), std::make_pair(&u2, "u2")};
  // std::vector<std::pair<ScrollingBuffer *, std::string>> measVecs{
  //     std::make_pair(&y0, "y0"),       std::make_pair(&y1, "y1"),
  //     std::make_pair(&y2, "y2"),       std::make_pair(&yhat0, "yhat0"),
  //     std::make_pair(&yhat1, "yhat1"), std::make_pair(&yhat2, "yhat2"),
  //     std::make_pair(&yref0, "yref0"), std::make_pair(&yref1, "yref1"),
  //     std::make_pair(&yref2, "yref2")};

  // template matching interactions
  // ImVector<ImVec2> points;
  // ImVec2 canvas_p0, canvas_p1, canvas_sz;
  // ImVec2 rectStart, rectEnd;
  // ImVec2 scrolling{0.0f, 0.0f};
  // bool opt_enable_grid = true;
  // bool opt_enable_rect = false;
  // bool opt_enable_context_menu = true;
  // bool adding_line = false;
  // bool addingRect = false;

  // for showing raw/processed frames
  // std::vector<cv::Mat> procFrames;
  // std::vector<int> procWidths, procHeights;

public:
  GUI(ImCap *imCap, ImProc *imProc, Pump *pump, Supervisor *sv);
  ~GUI();
  void startGUIThread();
  void imguiConfig();
  void imguiStyle();
  void render();
  void renderMenu();

  // void showImProcSetup();
  // void showImProc();
  // void showCtrlSetup();
  // void showCtrl();
  // void showSysIDSetup();
  // void showSysID();

  std::thread guiThread;
};
