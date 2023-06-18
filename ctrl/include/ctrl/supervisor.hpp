#pragma once

#include "SupervisoryController.h" // Model header file
#include "improc/improc.hpp"
#include "pump/pump.hpp"
#include "util/util.hpp"
#include <numeric>
#include <pybind11/embed.h>

class Supervisor {
  std::atomic<bool> startedCtrl{false}, startedSysIDFlag{false};
  std::thread ctrlThread, sysIDThread;

  // Called within thread context
  void start();

public:
  ordered_value conf;
  bool simModeActive;
  std::string dataPath, confPath;

  Pump *pump;
  ImProc *imProc;

  std::vector<double> y;

  SupervisoryController *sup;
  SupervisoryController::ExtU supIn;
  SupervisoryController::ExtY supOut;

  // defined in SupervisoryControler_data.cpp in SupervisoryController::P SupervisoryController::rtP
  constexpr static event_bus nullEv{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 0.0, 0.0, 0.0};
  event_bus currEv_;
  std::mutex currEvMtx;
  QueueFPS<event_bus> *evQueue_;
  real_T y_max[3], y_o[3], u_o[3], y_range[3];
  boolean_T inTransRegion;
  bool measAvail_;
  bool trueMeasAvail[3];
  time_point<steady_clock> initTime{steady_clock::now()};
  time_point<steady_clock> prevCtrlTime = initTime;

  QueueFPS<int> *ctrlDataQueuePtr;

  Supervisor(ImProc *imProc, Pump *pump);
  ~Supervisor();

  void startThread();
  void stopThread();

  void addEvent(event_bus e) { evQueue_->push_back(e); }

  void setCurrEv(event_bus e) {
    std::lock_guard<std::mutex> lock(currEvMtx);
    currEv_ = e;
  }
  event_bus getCurrEv() {
    std::lock_guard<std::mutex> lock(currEvMtx);
    return currEv_;
  }

  std::string getDataPath() const { return dataPath; }
  std::string getConfPath() const { return confPath; }
  void clearCtrlDataQueue();
};
