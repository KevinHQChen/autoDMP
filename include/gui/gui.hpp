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
#include "util/util.hpp"

#include <cstdio>

class GUI {
  ordered_value conf;
  guiConfig guiConf;

  ImCap *imCap = nullptr;
  // ImProc *imProc = nullptr;
  // Supervisor *supervisor = nullptr;

  GLFWwindow *window;
  GLuint textureID;
  ImGuiWindowFlags imCapFlags =
      ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
  bool needToQuit{false};
  std::optional<std::pair<int, int>> newSize{};

  cv::Mat rawFrame;
  int rawWidth, rawHeight;
  std::vector<QueueFPS<cv::Mat> *> tempResultsQueues, procFramesQueues;
  std::vector<QueueFPS<cv::Point> *> procDataQueues;

public:
  GUI();
  ~GUI();
  void startGUIThread();
  ImGuiWrapperReturnType render();
  int imguiMain();

  void showRawImCap();
  // void showProcImCap(bool &startImCap);
  void updateTexture(const cv::Mat &img);

  std::thread guiThread;
};
