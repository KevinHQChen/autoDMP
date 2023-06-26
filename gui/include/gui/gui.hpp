#pragma once

#include "gui/windows.hpp"
#include "imgui_md_wrapper.h"
#include "immapp/immapp.h"
#include "implot/implot.h"

#include "ctrl/supervisor.hpp"
#include "imcap/imcap.hpp"
#include "improc/improc.hpp"
#include "pump/pump.hpp"
#include "util/util.hpp"

#include <cmath>
#include <cstdio>

#include <spdlog/async.h>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/sinks/basic_file_sink.h>

namespace py = pybind11;
using namespace py::literals;

// Custom sink for spdlog that stores messages in a deque
template <typename Mutex> class ImGuiConsoleSink : public sinks::base_sink<Mutex> {
protected:
  void sink_it_(const details::log_msg &msg) override {
    // Store message in deque
    memory_buf_t formatted;
    sinks::base_sink<Mutex>::formatter_->format(msg, formatted);
    logBuffer.push_back(fmt::to_string(formatted));
  }

  void flush_() override {
    // Flush could be implemented here
  }

public:
  std::deque<std::string> logBuffer;
};

// Convenient typedefs for the console sink
using ImGuiConsoleSink_mt = ImGuiConsoleSink<std::mutex>;
using ImGuiConsoleSink_st = ImGuiConsoleSink<details::null_mutex>;

// The ImGui console class
class ImGuiConsole {
public:
  ImGuiConsole() {
    // // Create a thread pool with 1 thread and a queue with room for 8192 log messages
    // spdlog::init_thread_pool(8192, 1);
    // auto thread_pool = spdlog::thread_pool();

    // // Create the sink and add it to the spdlog
    // consoleSink = std::make_shared<ImGuiConsoleSink_mt>();
    // auto async_logger = std::make_shared<spdlog::async_logger>(
    //     "async_logger", consoleSink, thread_pool, spdlog::async_overflow_policy::block);
    // spdlog::set_default_logger(async_logger);

    // Create the sink and add it to the spdlog
    consoleSink = std::make_shared<ImGuiConsoleSink_mt>();
    default_logger()->sinks().push_back(consoleSink);
    fileSink = std::make_shared<sinks::basic_file_sink_mt>("log.txt", true);
    default_logger()->sinks().push_back(fileSink);
  }

  void render() {
    // Create an ImGui window
    ImGui::Begin("Console");

    // Display all the log messages
    for (const auto &message : consoleSink->logBuffer)
      ImGui::TextUnformatted(message.c_str());

    // Auto scroll to the bottom
    if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY())
      ImGui::SetScrollHereY(1.0f);

    ImGui::End();
  }

private:
  std::shared_ptr<ImGuiConsoleSink_mt> consoleSink;
  std::shared_ptr<sinks::basic_file_sink_mt> fileSink;
};

class GUI {
  HelloImGui::RunnerParams runnerParams;
  ImmApp::AddOnsParams addOnsParams;

  ordered_value conf;
  guiConfig guiConf;

  ImCap *imCap_;
  ImProc *imProc_;
  Pump *pump_;
  Supervisor *sv_;

  std::shared_ptr<gui::PumpWindow> pumpWindow_;
  std::shared_ptr<gui::ImProcWindow> imProcWindow_;
  std::shared_ptr<gui::CtrlWindow> ctrlWindow_;
  std::shared_ptr<gui::PlotWindow> plotWindow_;
  ImGuiConsole &console_;

public:
  GUI(ImCap *imCap, ImProc *imProc, Pump *pump, Supervisor *sv, ImGuiConsole &consoleSink);
  ~GUI();
  void startThread();
  void imguiConfig();
  void imguiStyle();
  void render();
  void renderMenu();

  std::thread guiThread;
};
