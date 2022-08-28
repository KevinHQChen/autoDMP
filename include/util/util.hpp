#pragma once

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
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

#include <toml++/toml.h>
#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <spdlog/spdlog.h>


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
