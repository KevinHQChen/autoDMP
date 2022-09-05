#pragma once

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

#include "imguiwrap.dear.h"
#include "imguiwrap.helpers.h"
#include <stdio.h>

#include "improc/imcap.hpp"
#include "improc/improc.hpp"
#include "util/util.hpp"

#include <cstdio>

class GUI {
  ordered_value conf;
  guiConfig guiConf;

  ImCap *imCap = nullptr;
  ImProc *imProc = nullptr;
  // Supervisor *supervisor = nullptr;

  GLFWwindow *window;
  GLuint rawTextureID, preTextureID;
  std::vector<GLuint> procTextureIDs;
  ImGuiWindowFlags imCapFlags =
      ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
  std::optional<std::pair<int, int>> newSize{};
  bool needToQuit{false};

  // for template matching
  ChannelPose chanPose;
  cv::Mat chans[4], firstChan, currChan;
  cv::Rect firstChanBBox, firstRotChanBBox, templateBBox;
  // template matching interactions
  ImVector<ImVec2> points;
  ImVec2 rectStart, rectEnd;
  ImVec2 scrolling{0.0f, 0.0f};
  bool opt_enable_grid = true;
  bool opt_enable_rect = false;
  bool opt_enable_context_menu = true;
  bool adding_line = false;
  bool addingRect = false;

  // for showing raw/processed frames
  cv::Mat rawFrame, preFrame;
  int rawWidth, rawHeight, preWidth, preHeight;
  std::vector<cv::Mat> procFrames;
  std::vector<int> procWidths, procHeights;

  std::vector<QueueFPS<cv::Point> *> procDataQueues;

public:
  GUI();
  ~GUI();
  void startGUIThread();
  ImGuiWrapperReturnType render();
  int imguiMain();

  void showRawImCap();
  void showTmplMatchSetup();
  void showImProc();
  void updateTexture(const cv::Mat &img, GLuint &textureID);

  std::thread guiThread;
};
