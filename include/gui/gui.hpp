#pragma once

#include "util/util.hpp"
#include "imgui.h"
#include "imgui_impl_opengl3.h"
#include "imgui_impl_sdl.h"
// alternative to SDL, lower level - more direct access to OpenGL
// #include "imgui_impl_glfw.h"
#include <GL/glew.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <GLFW/glfw3.h>

// #include "../include/downloader.h"

#include <stdio.h>

class GUIRenderer
{
public:
GUIRenderer();
~GUIRenderer(void);
int InitGUI();
void Render();
void ShowImage();
void UpdateTexture();

private:
SDL_Window *window_;
GLuint texture_id_;
int image_width_, image_height_;
// Downloader downloader_;
cv::Mat resized_image_;
cv::Mat thresholded_image_;
int threshold_;
};
