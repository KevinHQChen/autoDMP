#pragma once

#include "SupervisoryController.h" // Model header file
#include "improc/improc.hpp"
#include "pump/pump.hpp"
#include "util/util.hpp"
#include <numeric>
#include <pybind11/embed.h>

void *getSupervisorParam(rtwCAPI_ModelMappingInfo *mmi, uint_T paramIdx);

class Supervisor {
  std::atomic<bool> startedCtrl{false}, startedSysIDFlag{false};
  std::thread ctrlThread, sysIDThread;

  ordered_value conf;
  std::string dataPath, confPath;
  bool simModeActive;

  std::vector<double> y;

  event_bus currEv_;
  std::mutex currEvMtx;

  time_point<steady_clock> initTime{steady_clock::now()};
  time_point<steady_clock> prevCtrlTime = initTime;

  // Called within thread context
  void start();

public:
  Supervisor(ImProc *imProc, Pump *pump);
  ~Supervisor();

  Pump *pump;
  ImProc *imProc;

  SupervisoryController *sup;
  SupervisoryController::ExtU supIn;
  SupervisoryController::ExtY supOut;
  SupervisoryController::RT_MODEL *rtM;

  // defined in SupervisoryController_data.cpp in SupervisoryController::P
  // SupervisoryController::rtP
  event_bus nullEv;

  void startThread();
  void stopThread();

  QueueFPS<event_bus> *evQueue_;
  void addEvent(event_bus e) { evQueue_->push_back(e); }
  void setCurrEv(event_bus e) {
    std::lock_guard<std::mutex> lock(currEvMtx);
    currEv_ = e;
  }
  event_bus getCurrEv() {
    std::lock_guard<std::mutex> lock(currEvMtx);
    return currEv_;
  }

  QueueFPS<int> *ctrlDataQueuePtr;
  void clearCtrlDataQueue();

  std::string getDataPath() const { return dataPath; }
  std::string getConfPath() const { return confPath; }
};
