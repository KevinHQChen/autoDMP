#pragma once

#include "SupervisoryController.h" // Model header file
#include "improc/improc.hpp"
#include "pump/pump.hpp"
#include "util/util.hpp"
#include <pybind11/embed.h>

struct Event {
  int srcState, destState;
  Eigen::Vector3d destPos; // as % of channel length
  Eigen::Vector3d vel;     // as px/second

  Event(int srcState, int destState, Eigen::Vector3d destPos, Eigen::Vector3d vel)
      : srcState(srcState), destState(destState), destPos(destPos), vel(vel) {}
};

struct StateData;
class State; // forward declaration

class SupervisoryControllerRecvData_event_busT : public RecvData_event_busT{
 public:
  void RecvData(event_bus* data, int32_T length, int32_T* status)
  {
    // Add receive data logic here
  }
};

class Supervisor {

  StateData *currStateData_;

  std::atomic<bool> startedCtrl{false}, startedSysIDFlag{false};
  std::thread ctrlThread, sysIDThread;

  // Called within thread context
  void start();
  void startSysID();

public:
  ordered_value conf;
  bool simModeActive;
  std::string dataPath, confPath;

  Pump *pump;
  ImProc *imProc;
  SupervisoryController *sup;
  SupervisoryControllerRecvData_event_busT nextEventRecvData_arg;
  SupervisoryController::ExtU *supIn;
  SupervisoryController::ExtY *supOut;
  State *currState_ = nullptr;
  Event *currEvent_ = nullptr;
  QueueFPS<Event *> *eventQueue_;
  event_bus *currEv_ = nullptr;
  QueueFPS<event_bus *> *evQueue_;
  QueueFPS<int> *ctrlDataQueuePtr;

  Supervisor(ImProc *imProc, Pump *pump);
  ~Supervisor();

  void startThread();
  void stopThread();
  void startSysIDThread(Eigen::Vector3d uref, bool *selChs, std::vector<float> minVals,
                        std::vector<float> maxVals, Eigen::MatrixXd &data);
  void stopSysIDThread();

  bool measAvail();
  bool updateInputs();

  void addEvent(int srcState, int destState, Eigen::Vector3d pos, Eigen::Vector3d vel);

  // tmpl methods must be defined in headers (https://stackoverflow.com/a/10632266)
  template <typename T> void updateState() {
    delete currState_;
    currState_ = new T(this);
  }

  template <typename T> void updateState(Eigen::Vector3d uref_) {
    delete currState_;
    currState_ = new T(this, uref_);
  }

  template <typename T> void updateState(Eigen::Vector3d uref_, int prevState) {
    delete currState_;
    currState_ = new T(this, prevState, uref_);
  }

  template <typename T>
  void updateState(Eigen::Vector3d uref, bool *selChs, std::vector<float> minVals,
                   std::vector<float> maxVals, Eigen::MatrixXd &data) {
    delete currState_;
    currState_ = new T(this, uref, selChs, minVals, maxVals, data);
  }

  std::string getDataPath() const { return dataPath; }
  std::string getConfPath() const { return confPath; }
  StateData getCurrStateData();
  void clearCtrlDataQueue();
};
