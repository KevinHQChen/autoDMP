#pragma once

#include "util/util.hpp"
#include <GL/glew.h>

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

    return *this;
  }
};

inline void GUIFrame::updateTexture() {
  // mimic opencv grayscale image in opengl by making each color channel the same
  if (mat.empty()) {
    error("Image is empty");
    return;
  }

  // get image data type
  GLenum imgDataType;
  if (mat.type() == CV_8UC1) {
    imgDataType = GL_UNSIGNED_BYTE;
    // info("Image is CV_8UC1");
  } else if (mat.type() == CV_16UC1) {
    imgDataType = GL_UNSIGNED_SHORT;
    // info("Image is CV_16UC1");
  } else {
    error("Unsupported image data type");
    return;
  }

  cv::Mat tmp;
  cv::merge(std::vector<cv::Mat>{mat, mat, mat}, tmp);

  // https://github.com/ocornut/imgui/issues/4628
  if (texture == 0) {
    // update texture
    // glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    // create opengl texture identifier
    glGenTextures(1, &texture);
  }
  glBindTexture(GL_TEXTURE_2D, texture);
  // setup filtering parameters for display
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  // upload pixels into texture
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glPixelStorei(
      GL_UNPACK_ALIGNMENT,
      1); // https://stackoverflow.com/a/45486871
          // (https://www.khronos.org/opengl/wiki/Common_Mistakes#Texture_upload_and_pixel_reads)
  glTexImage2D(GL_TEXTURE_2D, // texture type
               0,             // pyramid level (for mip-mapping), 0 is the top level
               GL_RGB,        // internal color format to convert to
               // (https://www.khronos.org/opengl/wiki/Image_Format)
               tmp.cols, tmp.rows, // image width, height
               0,                  // border width in pixels (can be 1 or 0)
               GL_BGR,             // input image format
               imgDataType, // image data type (https://www.khronos.org/opengl/wiki/OpenGL_Type)
               tmp.ptr());  // pointer to data
  // glGenerateMipmap(GL_TEXTURE_2D);
}
