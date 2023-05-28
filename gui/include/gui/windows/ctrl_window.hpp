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

  event_bus currEv_;
  int srcState, destState, moveTime, holdTime;
  int targetPos[NUM_CHANS];
  bool chs0[NUM_CHANS] = {true, false, false}, chs1[NUM_CHANS] = {false, true, true},
       chs2[NUM_CHANS] = {true, false, true};
  float excitationAmp, paramLowerBound, covModification;

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

  void renderAddEventDialog();
  void renderEventQueueContents();
  void renderSupervisorStatus();

public:
  bool ctrlSetupVisible_{false}, ctrlVisible_{false};

  CtrlWindow(Supervisor *sv);
  ~CtrlWindow();
  void render() override;
};

} // namespace gui
