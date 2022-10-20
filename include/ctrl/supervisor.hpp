#pragma once

#include "improc/improc.hpp"
#include "pump/pump.hpp"
#include "util/util.hpp"

struct Event {
  int srcState, destState;
  Eigen::Vector3d destPos; // as % of channel length
  Eigen::Vector3d vel;     // as px/second

  Event(int srcState, int destState, Eigen::Vector3d destPos, Eigen::Vector3d vel)
      : srcState(srcState), destState(destState), destPos(destPos), vel(vel) {}
};

struct StateData;
class State; // forward declaration

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

  Pump *pump = nullptr;
  ImProc *imProc = nullptr;
  State *currState_ = nullptr;
  Event *currEvent_ = nullptr;
  QueueFPS<Event *> *eventQueue_;
  QueueFPS<int> *ctrlDataQueuePtr;

  bool gotPrbs = false;
  Eigen::MatrixXd prbs, prbs_, scaledPrbs;
  float prbsScale[3] = {1, 1, 1};
  const char *prbsRows[3] = {"row0", "row1", "row2"};
  float prbsUrefArr[3] = {60, 40, 60}; // default to state0 uref
  Eigen::Vector3d prbsUref = Eigen::Vector3d::Zero();

  Supervisor(ImProc *imProc, Pump *pump);
  ~Supervisor();

  void startThread();
  void stopThread();
  void startSysIDThread();
  void stopSysIDThread();
  bool started();
  bool startedSysID();

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


  std::string getDataPath() const { return dataPath; }
  std::string getConfPath() const { return confPath; }
  StateData getCurrStateData();
  void clearCtrlDataQueue();
};
