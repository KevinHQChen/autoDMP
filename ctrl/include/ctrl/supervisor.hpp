#pragma once

#include "SupervisoryController.h" // Model header file
#include "improc/improc.hpp"
#include "pump/pump.hpp"
#include "util/util.hpp"
#include <numeric>
#include <pybind11/embed.h>

struct Event {
  int srcState, destState;
  Eigen::Vector3d destPos; // as % of channel length
  Eigen::Vector3d vel;     // as px/second

  Event(int srcState, int destState, Eigen::Vector3d destPos, Eigen::Vector3d vel)
      : srcState(srcState), destState(destState), destPos(destPos), vel(vel) {}
};

class State; // forward declaration

class Supervisor {
  std::atomic<bool> startedCtrl{false}, startedSysIDFlag{false};
  std::thread ctrlThread, sysIDThread;

  // Called within thread context
  void start();
  void startSysID();
  void updateMeas();

public:
  ordered_value conf;
  bool simModeActive;
  std::string dataPath, confPath;

  Pump *pump;
  ImProc *imProc;

  std::vector<Pose> p;

  SupervisoryController *sup;
  SupervisoryController::ExtU supIn;
  SupervisoryController::ExtY supOut;

  // defined in SupervisoryControler_data.cpp in SupervisoryController::P SupervisoryController::rtP
  constexpr static event_bus nullEv{0.0,
                                    0.0,

                                    {0.0, 0.0, 0.0},
                                    0.0,
                                    0.0,

                                    {true, false, false},

                                    {true, false, false}};
  event_bus currEv_ = nullEv;
  QueueFPS<event_bus> *evQueue_;
  real_T y[3], y_max[3], y_o[3], u_o[3], y_range[3];
  boolean_T inTransRegion;
  bool measAvail_;
  bool trueMeasAvail[3];
  time_point<steady_clock> initTime{steady_clock::now()};
  time_point<steady_clock> prevCtrlTime = initTime;

  State *currState_ = nullptr;
  Event *currEvent_ = nullptr;
  QueueFPS<int> *ctrlDataQueuePtr;

  Supervisor(ImProc *imProc, Pump *pump);
  ~Supervisor();

  void startThread();
  void stopThread();
  void startSysIDThread(Eigen::Vector3d uref, bool *selChs, std::vector<float> minVals,
                        std::vector<float> maxVals, Eigen::MatrixXd &data);
  void stopSysIDThread();

  void addEvent(event_bus e) { evQueue_->push_back(e); }

  std::string getDataPath() const { return dataPath; }
  std::string getConfPath() const { return confPath; }
  void clearCtrlDataQueue();

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
};
