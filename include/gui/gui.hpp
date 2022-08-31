#pragma once

#include "util/util.hpp"

#include "imgui.h"
#include "imgui_impl_opengl3.h"
// glfw (more direct access to OpenGL) is a lower level alternative to SDL
#include "imgui_impl_glfw.h"
// SDL is more suited for 2D game dev
// #include "imgui_impl_sdl.h"

// renderers
#include <GL/glew.h>

// platforms
#include <GLFW/glfw3.h>
// #include <SDL2/SDL.h>
// #include <SDL2/SDL_opengl.h>

#include <stdio.h>
#include "imguiwrap.dear.h"
#include "imguiwrap.helpers.h"

#if 0
class GUIRenderer {
public:
  GUIRenderer();
  ~GUIRenderer(void);
  int InitGUI();
  void Render();
  void ShowImage();
  void UpdateTexture();

private:
  ImVec4 clear_color = ImVec4(0.45f, 0.56f, 0.67f, 1.00f);
  GLFWwindow *window_;
  GLuint texture_id_;
  int image_width_, image_height_;
  cv::Mat resized_image_;
  cv::Mat thresholded_image_;
  int threshold_;
};
#endif
