#pragma once

#include "imgui.h"
#include "implot.h"

// renderers
#include "imgui_impl_opengl3.h"
#include <GL/glew.h>

// platforms
// glfw (more direct access to OpenGL) is a lower level alternative to SDL
#include "imgui_impl_glfw.h"
#include <GLFW/glfw3.h>
// SDL is more suited for 2D game dev
// #include "imgui_impl_sdl.h"
// #include <SDL2/SDL.h>
// #include <SDL2/SDL_opengl.h>

// #include "imguiwrap.dear.h"
// #include "imguiwrap.helpers.h"
#include <stdio.h>

#include "ctrl/supervisor.hpp"
#include "gui/guiframe.hpp"
#include "improc/imcap.hpp"
#include "improc/improc.hpp"
#include "pump/pump.hpp"
#include "util/util.hpp"

#include <cstdio>

// utility structure for realtime plot
struct ScrollingBuffer {
  int MaxSize;
  int Offset;
  ImVector<ImVec2> Data;
  ScrollingBuffer(int max_size = 2000) {
    MaxSize = max_size;
    Offset = 0;
    Data.reserve(MaxSize);
  }
  void AddPoint(float x, float y) {
    if (Data.size() < MaxSize)
      Data.push_back(ImVec2(x, y));
    else {
      Data[Offset] = ImVec2(x, y);
      Offset = (Offset + 1) % MaxSize;
    }
  }
  void Erase() {
    if (Data.size() > 0) {
      Data.shrink(0);
      Offset = 0;
    }
  }
};

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
};

void setWindowFullscreen();

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
  ordered_value conf;
  guiConfig guiConf;

  ImCap *imCap = nullptr;
  ImProc *imProc = nullptr;
  Pump *pump = nullptr;
  Supervisor *sv = nullptr;

  GLFWwindow *window;
  std::optional<std::pair<int, int>> newSize{};
  std::vector<GLuint> procTextureIDs;
  ImGuiDockNodeFlags dockNodeFlags = ImGuiDockNodeFlags_PassthruCentralNode;
  // We are using the ImGuiWindowFlags_NoDocking flag to make the parent window not dockable into,
  // because it would be confusing to have two docking targets within each others.
  // When using ImGuiDockNodeFlags_PassthruCentralNode, DockSpace() will render our background
  // and handle the pass-thru hole, so we ask Begin() to not render a background.
  ImGuiWindowFlags dockSpaceFlags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking |
                                    ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse |
                                    ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                                    ImGuiWindowFlags_NoBringToFrontOnFocus |
                                    ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoBackground;
  ImGuiWindowFlags imCapFlags = 0;
  // ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
  ImGuiWindowFlags imProcSetupFlags = 0;
  bool needToQuit{false};

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
  GUIFrame rawFrame, preFrame;
  GUIFrame procGUIFrames[3];
  std::array<GUIFrame, NUM_TEMPLATES> tmplGUIFrames;
  std::vector<cv::Mat> procFrames;
  std::vector<int> procWidths, procHeights;

public:
  GUI();
  ~GUI();
  void startGUIThread();
  std::optional<int> render();
  int imguiMain();

  void showRawImCap();
  void showImProc();
  void showImProcSetup();
  void showCtrl();
  void showCtrlSetup();
  void showSysIDSetup();
  void showSysID();

  void contextMenu(bool enable);

  std::thread guiThread;
};
