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

  std::shared_ptr<gui::PumpWindow> pumpWindow_;
  std::shared_ptr<gui::ImProcWindow> imProcWindow_;
  std::shared_ptr<gui::CtrlWindow> ctrlWindow_;
  std::shared_ptr<gui::PlotWindow> plotWindow_;

public:
  GUI(ImCap *imCap, ImProc *imProc, Pump *pump, Supervisor *sv);
  ~GUI();
  void startGUIThread();
  void imguiConfig();
  void imguiStyle();
  void render();
  void renderMenu();

  std::thread guiThread;
};
