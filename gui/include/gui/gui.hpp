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
#include <spdlog/sinks/stdout_color_sinks.h>

namespace py = pybind11;
using namespace py::literals;

// Custom sink for spdlog that stores messages in a deque
class ImGuiConsoleSinkMt : public sinks::base_sink<std::mutex> {
public:
  ImGuiConsoleSinkMt() {
    // Set a default pattern formatter
    set_formatter(std::make_unique<spdlog::pattern_formatter>());
  }

  std::deque<std::string> logBuffer;

protected:
  void sink_it_(const details::log_msg &msg) override {
    // Store formatted message in deque
    spdlog::memory_buf_t formatted;
    formatter_->format(msg, formatted);
    logBuffer.push_back(fmt::to_string(formatted));
  }

  void flush_() override {
    // Flush could be implemented here
  }
};

// The ImGui console class
class ImGuiConsole {
public:
  ImGuiConsole() {
    consoleSink = std::make_shared<ImGuiConsoleSinkMt>();
    fileSink = std::make_shared<sinks::basic_file_sink_mt>("log.txt", true);
    stdOutSink = std::make_shared<sinks::stdout_color_sink_mt>();
  }

  std::shared_ptr<logger> getLogger(const std::string &name) {
    auto logger_ =
        std::make_shared<logger>(name, sinks_init_list{consoleSink, fileSink, stdOutSink});
    register_logger(logger_);
    return logger_;
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
  std::shared_ptr<ImGuiConsoleSinkMt> consoleSink;
  std::shared_ptr<sinks::basic_file_sink_mt> fileSink;
  std::shared_ptr<sinks::stdout_color_sink_mt> stdOutSink;
};

class GUI {
public:
  GUI(ImCap *imCap, ImProc *imProc, Pump *pump, Supervisor *sv, ImGuiConsole &consoleSink,
      std::shared_ptr<logger> log);
  ~GUI();
  void startThread();
  void imguiConfig();
  void imguiStyle();
  void render();
  void renderMenu();

  std::thread guiThread;

private:
  std::shared_ptr<logger> lg;
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
};
