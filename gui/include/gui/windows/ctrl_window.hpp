#pragma once

#include "gui/components/button.hpp"
#include "gui/components/implot_helpers.hpp"
#include "gui/components/slider.hpp"
#include "window.hpp"

#include "ctrl/supervisor.hpp"

namespace gui {

class CtrlWindow : public Window {
  Supervisor *sv_;

  // displaying/modifying events, states
  ImGuiTableFlags tableFlags = ImGuiTableFlags_BordersV | ImGuiTableFlags_BordersOuterH |
                               ImGuiTableFlags_Resizable | ImGuiTableFlags_RowBg |
                               ImGuiTableFlags_NoBordersInBody;

  float guiTime{0.0f}, history{30.0f};
  ScrollingBuffer u0, u1, u2;
  ScrollingBuffer y0, y1, y2, yhat0, yhat1, yhat2, yref0, yref1, yref2;
  ScrollingBuffer param0, param1, param2, param3, param4, param5, param6, param7, param8;
  std::vector<std::pair<ScrollingBuffer *, std::string>> ctrlVecs{
      std::make_pair(&u0, "u0"), std::make_pair(&u1, "u1"), std::make_pair(&u2, "u2")};
  std::vector<std::pair<ScrollingBuffer *, std::string>> measVecs{
      std::make_pair(&y0, "y0"),       std::make_pair(&y1, "y1"),
      std::make_pair(&y2, "y2"),       std::make_pair(&yhat0, "yhat0"),
      std::make_pair(&yhat1, "yhat1"), std::make_pair(&yhat2, "yhat2"),
      std::make_pair(&yref0, "yref0"), std::make_pair(&yref1, "yref1"),
      std::make_pair(&yref2, "yref2")};
  std::vector<std::pair<ScrollingBuffer *, std::string>> paramVecs{
      std::make_pair(&param0, "param0"), std::make_pair(&param1, "param1"),
      std::make_pair(&param2, "param2"), std::make_pair(&param3, "param3"),
      std::make_pair(&param4, "param4"), std::make_pair(&param5, "param5"),
      std::make_pair(&param6, "param6"), std::make_pair(&param7, "param7"),
      std::make_pair(&param8, "param8")};
  int srcState, destState, moveTime, holdTime;
  int targetPos[NUM_CHANS];
  bool chs0[NUM_CHANS] = {true, false, false}, chs1[NUM_CHANS] = {false, true, true},
       chs2[NUM_CHANS] = {true, false, true};
  float excitationAmp;

  int openAction = -1;
  int dropletLength = 0;
  int numDroplets = 0;

  event_bus getEvent(int srcState, int destState, std::array<int, NUM_CHANS> targetPos,
                     int moveTime, int holdTime) {
    event_bus e;
    e.srcState = srcState;
    e.destState = destState;
    for (int i = 0; i < NUM_CHANS; ++i)
      e.destPos[i] = targetPos[i] / 100.0;
    e.moveTime = moveTime;
    e.holdTime = holdTime;

    switch (srcState) {
    case 0: // [1 0 0]
      std::memcpy(e.chs, chs0, sizeof(chs0));
      break;
    case 1: // [0 1 1]
      std::memcpy(e.chs, chs1, sizeof(chs1));
      break;
    case 2: // [1 0 1]
      std::memcpy(e.chs, chs2, sizeof(chs2));
      break;
    }
    switch (destState) {
    case 0: // [1 0 0]
      std::memcpy(e.nextChs, chs0, sizeof(chs0));
      break;
    case 1: // [0 1 1]
      std::memcpy(e.nextChs, chs1, sizeof(chs1));
      break;
    case 2: // [1 0 1]
      std::memcpy(e.nextChs, chs2, sizeof(chs2));
      break;
    }
    return e;
  }

public:
  bool ctrlSetupVisible_{false}, ctrlVisible_{false}, pauseCtrlDataViz{false};

  CtrlWindow(Supervisor *sv);
  ~CtrlWindow();
  void render() override;
};

} // namespace gui
