#pragma once

#include "improc/improc.hpp"
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
  // QueueFPS<Eigen::Matrix<int16_t, 3, 1>> *ctrlDataQueuePtr;

  std::atomic<bool> startedCtrl{false};
  std::thread ctrlThread;

  // Called within ctrlThread context
  void start();

public:
  ordered_value conf;
  std::string dataPath, confPath;

  ImProc *imProc = nullptr;
  State *currState_ = nullptr;
  Event *currEvent_ = nullptr;
  QueueFPS<Event *> *eventQueue_;

  Supervisor(ImProc *imProc);
  ~Supervisor();

  void startThread();
  void stopThread();
  bool started();

  void addEvent(int srcState, int destState, Eigen::Vector3d pos, Eigen::Vector3d vel);

  // tmpl methods must be defined in headers (https://stackoverflow.com/a/10632266)
  template <typename T> void updateState() {
    delete currState_;
    currState_ = new T(this);
  }

  std::string getDataPath() const { return dataPath; }
  std::string getConfPath() const { return confPath; }
  StateData getCurrStateData();
  void clearCtrlDataQueue();
};
