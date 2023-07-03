#pragma once

#include "SupervisoryController.h" // Model header file
#include "improc/improc.hpp"
#include "pump/pump.hpp"
#include "util/util.hpp"
#include <numeric>
#include <pybind11/embed.h>

class Supervisor {
public:
  Supervisor(ImProc *imProc, Pump *pump, std::shared_ptr<logger> log);
  ~Supervisor();

  Pump *pump;
  ImProc *imProc;

  SupervisoryController *sup;
  SupervisoryController::ExtU supIn;
  SupervisoryController::ExtY supOut;
  SupervisoryController::RT_MODEL *rtM;

  int no;

  event_bus nullEv;

  void startThread();
  void stopThread();
  bool started();

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

private:
  std::shared_ptr<logger> lg;
  std::atomic<bool> startedCtrl{false}, startedSysIDFlag{false};
  std::thread ctrlThread, sysIDThread;

  ordered_value conf;
  std::string dataPath, confPath;
  bool simModeActive;

  std::vector<double> y;
  std::vector<unsigned char> yo;

  event_bus currEv_;
  std::mutex ctrlMtx, currEvMtx;

  // Called within thread context
  void start();

  /**
   * Function: getSupervisorParam
   * ----------------------------
   * Retrieves the address of indexed variable block parameters defined in
   * `$GITROOT/external/SupervisoryController/scripts/SupervisoryController_ert_rtw/SupervisoryController_data.cpp`.
   *
   * @param mmi: A pointer to the model mapping information structure.
   * @param paramIdx: The index of the parameter whose address is to be retrieved.
   *
   * @return: A void pointer to the address of the specified parameter. This must be cast to the
   * appropriate type. If the model parameters are not available or the parameter address is not
   * available, the function returns nullptr.
   */
  void *getSupervisorParam(rtwCAPI_ModelMappingInfo *mmi, uint_T paramIdx);
};
