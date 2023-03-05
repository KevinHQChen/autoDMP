#pragma once

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

  std::shared_ptr<Pump> pump;
  std::shared_ptr<ImProc> imProc;
  State *currState_ = nullptr;
  Event *currEvent_ = nullptr;
  QueueFPS<Event *> *eventQueue_;
  QueueFPS<int> *ctrlDataQueuePtr;

  Supervisor(std::shared_ptr<ImProc> imProc, std::shared_ptr<Pump> pump);
  ~Supervisor();

  void startThread();
  void stopThread();
  void startSysIDThread(Eigen::Vector3d uref, bool *selChs, std::vector<float> minVals,
                        std::vector<float> maxVals, Eigen::MatrixXd &data);
  void stopSysIDThread();

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