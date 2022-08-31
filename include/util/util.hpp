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

// usage (print variable type in human-readable format): info("var type: {}",
// type_name<decltype(var)>());
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

struct config {
  CLI::App cli{fmt::format("{} version {}", "project name", "0.0.1")};
  // app = { fmt::format("{} version {}", autoDMP::cmake::project_name,
  // autoDMP::cmake::project_version) };

  CLI::App *cameraCli;

  bool paramsValid;     // 0:Invalid,   1:Valid
  int videoSource;      // 0:Webcam,    1:From file,       2:Andor
  std::string vidSrcFn; // used only if videoSource = 1
  int templateSource;   // 0:From file, 1:From videoSource
  int chanSource;       // 0:From file, 1:From videoSource
  // template matching threshold percentage (8-bit, i.e. pixel intensity/255
  // in
  // %) (0.5 is good for our current offline videos, 0.75 is good for live
  // videos)
  double tmThres; // replace this with repl later
  int saveRaw;    // 0:No,        1:Yes
  double rawFPS;  // used only if saveRaw = 1
  int saveProc;   // 0:No,        1:Yes
  double procFPS; // used only if saveProc = 1
  int imProc;     // 0:No,        1:Yes
  int ctrl;       // 0:No,        1:Yes
  int useGPU;     // 0:No,        1:Yes
  int numChannels;
  int dryRun;

  config();

  int parse(int argc, const char **argv);
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
