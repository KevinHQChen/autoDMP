#pragma once

#include <array>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgcodecs.hpp" // image formats to save to
#include "opencv2/imgproc.hpp"
#include "opencv2/opencv.hpp"

#include <Eigen/Dense>

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <spdlog/spdlog.h>
#include <toml.hpp>
#include <tsl/ordered_map.h>

// #include <string_view>
#include <boost/pfr.hpp>
#include <type_traits>

using namespace spdlog;
using namespace std::chrono;

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

template <typename T> struct isArray : std::false_type {};

template <typename T, std::size_t N> struct isArray<std::array<T, N>> : std::true_type {};

template <typename T> std::vector<std::variant<int, double>> toVector(const T &s) {
  std::vector<std::variant<int, double>> v;
  boost::pfr::for_each_field(s, [&](const auto &field) {
    if constexpr (isArray<std::decay_t<decltype(field)>>::value)
      for (const auto &item : field)
        v.push_back(item);
    else
      v.push_back(field);
  });
  return v;
}

template <typename T> T toStruct(const std::vector<std::variant<int, double>> &v) {
  T s;
  size_t i = 0;
  boost::pfr::for_each_field(s, [&](auto &field) {
    if constexpr (isArray<std::decay_t<decltype(field)>>::value) {
      for (auto &item : field) {
        if constexpr (std::is_same_v<decltype(item), int &>)
          item = std::get<int>(v[i]);
        else if constexpr (std::is_same_v<decltype(item), double &>)
          item = std::get<double>(v[i]);
        ++i;
      }
    } else {
      if constexpr (std::is_same_v<decltype(field), int &>)
        field = std::get<int>(v[i]);
      else if constexpr (std::is_same_v<decltype(field), double &>)
        field = std::get<double>(v[i]);
      ++i;
    }
  });
  return s;
}

class Timer {
public:
  Timer(const std::string &label) : start(high_resolution_clock::now()), label_(label) {}

  ~Timer() {
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    info("{} duration: {} ms", label_, duration.count());
  }

private:
  high_resolution_clock::time_point start;
  std::string label_;
};

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
    showDebug = toml::find<bool>(v, "showDebug");
  }
};

namespace Config {
static ordered_value conf = TOML11_PARSE_IN_ORDER("config/setup.toml");
static guiConfig guiConf = toml::find<guiConfig>(conf, "gui");
} // namespace Config

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

template <typename T> class SharedBuffer {
public:
  T get() {
    std::unique_lock<std::mutex> lock(mtx);
    cond_var.wait(lock);
    return item;
  }

  void set(const T &newItem) {
    std::lock_guard<std::mutex> lock(mtx);
    item = newItem;
    cond_var.notify_all();
  }

  void clear() {
    std::lock_guard<std::mutex> lock(mtx);
    item = T();
  }

private:
  T item;
  std::mutex mtx;
  std::condition_variable cond_var;
};

// create a generic type queue class template for storing frames
template <typename T> class QueueFPS : public std::deque<T> {
  std::string fileName;
  cv::TickMeter tm, tmSinceLastPush;
  std::mutex mutex;
  unsigned int counter;

public:
  // constructor:
  // chrono uses the concepts of timepoints and durations
  // now() provides a single timepoint (num of ticks since epoch) from a
  // specific clock this has no meaning by itself (unless you want ticks since
  // epoch, in which case use now().time_since_epoch().count()) a difference
  // between 2 timepoints returns a duration for which the number of ticks is
  // given by duration.count()
  QueueFPS(std::string fileName)
      : counter(0), fileName(fileName), out(fileName), startTime(std::chrono::steady_clock::now()) {
  }

  void push(const T &entry) {
    std::lock_guard<std::mutex> lockGuard(mutex);
    this->push_back(entry);
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
    this->pop_front();
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
      this->pop_front();
  }

  void clearFile() {
    std::lock_guard<std::mutex> lockGuard(mutex);
    out.fileStream.close();
    std::remove(fileName.c_str());
    out.fileStream.open(fileName);
    if (!out.fileStream.is_open()) {
      error("Could not open file");
    }
    out.fileStream << fileName << "\n";
    out.fileStream.flush(); // "commits" current buffer to file (in case close() is not called)
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
};
