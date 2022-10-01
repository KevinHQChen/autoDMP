#pragma once

#include "util/util.hpp"
#include "improc/improc.hpp"

struct Event {
  int srcState, destState;
  Eigen::Vector3d destPos; // as % of channel length
  Eigen::Vector3d vel; // as px/second

  Event(int srcState, int destState, Eigen::Vector3d destPos, Eigen::Vector3d vel) :
    srcState(srcState), destState(destState), destPos(destPos), vel(vel) {}
};

struct StateData;
class State; // forward declaration

class Supervisor {
  ordered_value conf;
  std::string dataPath;

  State *currState_ = nullptr;
  StateData *currStateData_;
  // QueueFPS<Eigen::Matrix<int16_t, 3, 1>> *ctrlDataQueuePtr;

  std::atomic<bool> startedCtrl{false};
  std::thread ctrlThread;

  // Called within ctrlThread context
  void start();

public:
  Event *currEvent_ = nullptr;
  QueueFPS<Event *> *eventQueue_;
  ImProc *imProc = nullptr;

  Supervisor(ImProc *imProc);
  ~Supervisor();

  void startThread();
  void stopThread();
  bool started();

  void addEvent(int srcState, int destState, Eigen::Vector3d pos, Eigen::Vector3d vel);

  template <typename T> void updateState();

  std::string getDataPath() const { return dataPath; }
  StateData getCurrStateData();
  void clearCtrlDataQueue();

};
