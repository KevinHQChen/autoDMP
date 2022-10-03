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
#include "util/util.hpp"

#include <cstdio>

// utility structure for realtime plot
struct ScrollingBuffer {
    int MaxSize;
    int Offset;
    ImVector<ImVec2> Data;
    ScrollingBuffer(int max_size = 2000) {
        MaxSize = max_size;
        Offset  = 0;
        Data.reserve(MaxSize);
    }
    void AddPoint(float x, float y) {
        if (Data.size() < MaxSize)
            Data.push_back(ImVec2(x,y));
        else {
            Data[Offset] = ImVec2(x,y);
            Offset =  (Offset + 1) % MaxSize;
        }
    }
    void Erase() {
        if (Data.size() > 0) {
            Data.shrink(0);
            Offset  = 0;
        }
    }
};

// utility structure for realtime plot
struct RollingBuffer {
    float Span;
    ImVector<ImVec2> Data;
    RollingBuffer() {
        Span = 10.0f;
        Data.reserve(2000);
    }
    void AddPoint(float x, float y) {
        float xmod = fmodf(x, Span);
        if (!Data.empty() && xmod < Data.back().x)
            Data.shrink(0);
        Data.push_back(ImVec2(xmod, y));
    }
};

void setWindowFullscreen();

class GUI {
  ordered_value conf;
  guiConfig guiConf;

  ImCap *imCap = nullptr;
  ImProc *imProc = nullptr;
  Supervisor *supervisor = nullptr;

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
  float guiTime{0.0f}, history{10.0f};
  ScrollingBuffer u0, u1, u2, y0, y1, y2, yref0, yref1, yref2;

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

  void contextMenu(bool enable);

  std::thread guiThread;
};
