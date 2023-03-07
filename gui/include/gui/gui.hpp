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

struct GUIEvent {
  int srcState = 0;
  int destState = 0;
  int pos[3] = {0, 0, 0};
  int vel[3] = {0, 0, 0};
  int *data[4];
  inline static const std::string props[4] = {"Src State", "Dest State", "Target Pos (ch0-2) [%]",
                                              "Target Vel (ch0-2) [px/s]"};
  inline static const int min[4] = {0, 0, 0, 0};
  inline static const int max[4] = {3, 3, 100, 20};

  GUIEvent() {
    int i = 0;
    data[i++] = &srcState;
    data[i++] = &destState;
    data[i++] = pos;
    data[i++] = vel;
  }

  GUIEvent(int srcState, int destState, Eigen::Vector3d pos, Eigen::Vector3d vel) {
    this->srcState = srcState;
    this->destState = destState;
    for (int i = 0; i < 3; i++) {
      this->pos[i] = pos(i);
      this->vel[i] = vel(i);
    }

    int i = 0;
    data[i++] = &srcState;
    data[i++] = &destState;
    data[i++] = this->pos;
    data[i++] = this->vel;
  }
};

inline void HelpMarker(const char *desc) {
  ImGui::TextDisabled("(?)");
  if (ImGui::IsItemHovered()) {
    ImGui::BeginTooltip();
    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
    ImGui::TextUnformatted(desc);
    ImGui::PopTextWrapPos();
    ImGui::EndTooltip();
  }
}

inline void displayVector3d(const char *vecName, Eigen::Vector3d vec) {
  ImGui::TableNextRow();
  ImGui::TableSetColumnIndex(0);
  ImGui::Text("%s", vecName);
  for (int i = 0; i < 3; ++i) {
    ImGui::TableSetColumnIndex(i + 1);
    ImGui::Text("%f", vec(i));
  }
}

inline void displayArray3b(const char *arrName, bool arr[3], const char *helpText = "") {
  ImGui::TableNextRow();
  ImGui::TableSetColumnIndex(0);
  ImGui::Text("%s", arrName);
  ImGui::SameLine();
  HelpMarker(helpText);
  for (int i = 0; i < 3; ++i) {
    ImGui::TableSetColumnIndex(i + 1);
    ImGui::Text("%d", arr[i]);
  }
}

class GUI {
  HelloImGui::RunnerParams runnerParams;
  ImmApp::AddOnsParams addOnsParams;

  ordered_value conf;
  guiConfig guiConf;

  std::shared_ptr<ImCap> imCap_;
  std::shared_ptr<ImProc> imProc_;
  std::shared_ptr<Pump> pump_;
  std::shared_ptr<Supervisor> sv_;

  std::shared_ptr<gui::SysIdWindow> sysIDWindow_;
  std::shared_ptr<gui::PumpWindow> pumpWindow_;
  std::shared_ptr<gui::ImProcWindow> imProcWindow_;

  // real time plotting
  ImPlotAxisFlags implotFlags = ImPlotAxisFlags_NoTickLabels;
  float guiTime{0.0f}, history{30.0f};
  ScrollingBuffer u0, u1, u2, du0, du1, du2;
  ScrollingBuffer y0, y1, y2, yref0, yref1, yref2;
  ScrollingBuffer dxhat0, dxhat1, dxhat2, z0, z1, z2;
  std::vector<std::pair<ScrollingBuffer *, std::string>> sysidCtrlVecs{
      std::make_pair(&u0, "u0"), std::make_pair(&u1, "u1"), std::make_pair(&u2, "u2")};
  std::vector<std::pair<ScrollingBuffer *, std::string>> sysidMeasVecs{
      std::make_pair(&y0, "y0"), std::make_pair(&y1, "y1"), std::make_pair(&y2, "y2")};
  std::vector<std::pair<ScrollingBuffer *, std::string>> ctrlVecs{
      std::make_pair(&u0, "u0"),   std::make_pair(&u1, "u1"),   std::make_pair(&u2, "u2"),
      std::make_pair(&du0, "du0"), std::make_pair(&du1, "du1"), std::make_pair(&du2, "du2")};
  std::vector<std::pair<ScrollingBuffer *, std::string>> measVecs{
      std::make_pair(&y0, "y0"),       std::make_pair(&y1, "y1"),
      std::make_pair(&y2, "y2"),       std::make_pair(&yref0, "yref0"),
      std::make_pair(&yref1, "yref1"), std::make_pair(&yref2, "yref2")};
  std::vector<std::pair<ScrollingBuffer *, std::string>> errorVecs{
      std::make_pair(&dxhat0, "dxhat0"), std::make_pair(&dxhat1, "dxhat1"),
      std::make_pair(&dxhat2, "dxhat2"), std::make_pair(&z0, "z0"),
      std::make_pair(&z1, "z1"),         std::make_pair(&z2, "z2")};
  void plotVector3d(const char *plotName, const char *xAx, const char *yAx, double yMin,
                    double yMax, std::vector<std::pair<ScrollingBuffer *, std::string>> &vecs);

  // displaying/modifying events, states
  ImGuiTableFlags tableFlags = ImGuiTableFlags_BordersV | ImGuiTableFlags_BordersOuterH |
                               ImGuiTableFlags_Resizable | ImGuiTableFlags_RowBg |
                               ImGuiTableFlags_NoBordersInBody;
  GUIEvent currEvent;
  std::deque<GUIEvent> guiEventQueue;
  int openAction = -1;
  int dropletLength = 0;

  // template matching interactions
  ImVector<ImVec2> points;
  ImVec2 canvas_p0, canvas_p1, canvas_sz;
  ImVec2 rectStart, rectEnd;
  ImVec2 scrolling{0.0f, 0.0f};
  bool opt_enable_grid = true;
  bool opt_enable_rect = false;
  bool opt_enable_context_menu = true;
  bool adding_line = false;
  bool addingRect = false;

  // for showing raw/processed frames
  std::vector<cv::Mat> procFrames;
  std::vector<int> procWidths, procHeights;

public:
  GUI(std::shared_ptr<ImCap> imCap, std::shared_ptr<ImProc> imProc, std::shared_ptr<Pump> pump,
      std::shared_ptr<Supervisor> sv);
  ~GUI();
  void startGUIThread();
  void imguiConfig();
  void imguiStyle();
  void render();
  void renderMenu();

  void showImProcSetup();
  void showImProc();
  void showCtrlSetup();
  void showCtrl();
  void showSysIDSetup();
  void showSysID();

  std::thread guiThread;
};
