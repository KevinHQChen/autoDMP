#pragma once

#include "imgui.h"

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

#include "improc/imcap.hpp"
#include "improc/improc.hpp"
#include "ctrl/supervisor.hpp"
#include "util/util.hpp"

#include <cstdio>

void setWindowFullscreen();

struct GUIFrame {
  cv::Mat mat;
  int width;
  int height;
  bool empty;
  GLuint texture{0};

  void updateTexture();

  GUIFrame &operator=(const cv::Mat &matInstance) {
    // grab most recent non-empty frame
    mat = matInstance;
    if (!mat.empty()) {
      empty = false;
      updateTexture();
      width = mat.cols;
      height = mat.rows;
      // or repeat previous frame if no new frame available
    } else if (width > 0 && height > 0)
      empty = false;
    else
      empty = true;

    return *this; // this is apparently inferred by compiler? still works without it
  }
};

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

  void contextMenu(bool enable);

  std::thread guiThread;
};
