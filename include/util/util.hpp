#pragma once

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgcodecs.hpp" // image formats to save to
#include "opencv2/imgproc.hpp"
#include "opencv2/opencv.hpp"

// #include "boost/program_options.hpp"
// #include "opencv2/core/cuda.hpp"
#include <Eigen/Dense>

// #include <CLI/CLI.hpp>
// #include <nfd.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <spdlog/spdlog.h>
#include <toml.hpp>
#include <tsl/ordered_map.h>

#include <string_view>

using namespace spdlog;
using namespace std::chrono;
using Vector1d = Eigen::Matrix<double, 1, 1>;
using Vector1ui = Eigen::Matrix<unsigned int, 1, 1>;
using Vector2ui = Eigen::Matrix<unsigned int, 2, 1>;

using ordered_value = toml::basic_value<toml::discard_comments, tsl::ordered_map, std::vector>;

#define TOML11_PARSE_IN_ORDER(...)                                                                 \
  toml::parse<toml::discard_comments, tsl::ordered_map>(__VA_ARGS__)

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
  bool startImProc;
  bool startImProcSetup;
  bool startPumpSetup;
  bool startCtrl;
  bool startCtrlSetup;
  bool startSysID;
  bool startSysIDSetup;
  bool pauseCtrlDataViz;
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
    startImProc = toml::find<bool>(v, "startImProc");
    startImProcSetup = toml::find<bool>(v, "startImProcSetup");
    startPumpSetup = toml::find<bool>(v, "startPumpSetup");
    startCtrl = toml::find<bool>(v, "startCtrl");
    startCtrlSetup = toml::find<bool>(v, "startCtrlSetup");
    startSysID = toml::find<bool>(v, "startSysID");
    startSysIDSetup = toml::find<bool>(v, "startSysIDSetup");
    pauseCtrlDataViz = toml::find<bool>(v, "pauseCtrlDataViz");
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
  st.fileStream.flush(); // "commits" current buffer to file (in case close() is not called)
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

// create a generic type queue class template for storing frames
template <typename T> class QueueFPS : public std::queue<T> {
public:
  // constructor:
  // chrono uses the concepts of timepoints and durations
  // now() provides a single timepoint (num of ticks since epoch) from a
  // specific clock this has no meaning by itself (unless you want ticks since
  // epoch, in which case use now().time_since_epoch().count()) a difference
  // between 2 timepoints returns a duration for which the number of ticks is
  // given by duration.count()
  QueueFPS(std::string fileName)
      : counter(0), out(fileName), startTime(std::chrono::steady_clock::now()) {}

  void push(const T &entry) {
    std::lock_guard<std::mutex> lockGuard(mutex);
    std::queue<T>::push(entry);
    out << std::fixed << std::setprecision(3);
    out << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() -
                                                                 startTime)
               .count()
        << ", ";
    counter += 1;
    if (counter == 1) {
      // Start counting from a second frame (warmup).
      tm.reset();
      tm.start();
    }
    tmSinceLastPush.reset();
    tmSinceLastPush.start();
  }

  T get() {
    std::lock_guard<std::mutex> lockGuard(mutex);
    T entry = this->front();
    this->pop();
    return entry;
  }

  float getFPS() {
    // std::lock_guard<std::mutex> lockGuard(mutex);
    tm.stop();
    double fps = counter / tm.getTimeSec();
    tm.start();
    return static_cast<float>(fps);
  }

  int getTimeSinceLastPush() {
    std::lock_guard<std::mutex> lockGuard(mutex);
    tmSinceLastPush.stop();
    double time = tmSinceLastPush.getTimeMilli();
    tmSinceLastPush.start();
    return (int)std::round(time);
  }

  void clear() {
    std::lock_guard<std::mutex> lockGuard(mutex);
    while (!this->empty())
      this->pop();
  }

  bool empty_() {
    std::lock_guard<std::mutex> lockGuard(mutex);
    return this->empty();
  }

  unsigned int counter_() {
    // std::lock_guard<std::mutex> lockGuard(mutex);
    return counter;
  }

  ffstream out;
  std::chrono::time_point<std::chrono::steady_clock> startTime;

private:
  cv::TickMeter tm, tmSinceLastPush;
  std::mutex mutex;
  unsigned int counter;
};
