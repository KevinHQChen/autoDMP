#pragma once

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgcodecs.hpp" // image formats to save to
#include "opencv2/imgproc.hpp"
#include "opencv2/opencv.hpp"

#include "boost/program_options.hpp"
#include "opencv2/core/cuda.hpp"
#include <Eigen/Dense>

#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <spdlog/spdlog.h>
#include <toml.hpp>
#include <tsl/ordered_map.h>

#include <string_view>

using namespace spdlog;
using namespace std::chrono;
using ordered_value = toml::basic_value<toml::discard_comments, tsl::ordered_map, std::vector>;

// print variable type in human-readable format
// usage (): info("var type: {}", type_name<decltype(var)>());
template <typename T> constexpr auto type_name() {
  std::string_view name, prefix, suffix;
#ifdef __clang__
  name = __PRETTY_FUNCTION__;
  prefix = "auto type_name() [T = ";
  suffix = "]";
#elif defined(__GNUC__)
  name = __PRETTY_FUNCTION__;
  prefix = "constexpr auto type_name() [with T = ";
  suffix = "]";
#elif defined(_MSC_VER)
  name = __FUNCSIG__;
  prefix = "auto __cdecl type_name<";
  suffix = ">(void)";
#endif
  name.remove_prefix(prefix.size());
  name.remove_suffix(suffix.size());
  return name;
}

struct guiConfig {
  std::string windowTitle;
  int width;
  int height;
  std::vector<float> clearColor;
  int enableVsync;
  bool keyboardNav;
  bool startDark;
  float fontSize;
  float scale;
  std::string fontPath;
  bool startImCap;
  bool showDebug;

  void from_toml(const ordered_value &v) {
    windowTitle = toml::find<std::string>(v, "windowTitle");
    width = toml::find<int>(v, "width");
    height = toml::find<int>(v, "height");
    clearColor = toml::find<std::vector<float>>(v, "clearColor");
    enableVsync = toml::find<int>(v, "enableVsync");
    keyboardNav = toml::find<bool>(v, "keyboardNav");
    startDark = toml::find<bool>(v, "startDark");
    fontSize = toml::find<float>(v, "fontSize");
    scale = toml::find<float>(v, "scale");
    fontPath = toml::find<std::string>(v, "fontPath");
    startImCap = toml::find<bool>(v, "startImCap");
    showDebug = toml::find<bool>(v, "showDebug");
  }
};

struct ffstream {
  std::ofstream fileStream;
  ffstream(std::string);
  ~ffstream();
};

template <class T> ffstream &operator<<(ffstream &st, T val) {
  st.fileStream << val;
  return st;
};

struct mstream {
  std::ofstream multiStream;
  mstream(std::string);
  ~mstream();
};

// Defining multiple output for the << operator
// 1 - to a file in the log folder
// 2 - to the console, as usual
template <class T> mstream &operator<<(mstream &st, T val) {
  st.multiStream << val;
  std::cout << val;
  return st;
};

void saveData(std::string fileName, Eigen::MatrixXd matrix);

Eigen::MatrixXd openData(std::string fileToOpen);
