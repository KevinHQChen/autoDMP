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
  int targetPos[2 * NUM_CHANS];
  float excitationAmp, paramLowerBound, covModification;

  int openAction = -1;
  int dropletLength = 0;
  int numDroplets = 0;

  event_bus getEvent(std::array<int, 2 * NUM_CHANS> targetPos, int moveTime, int holdTime) {
    event_bus e;
    for (int i = 0; i < 2 * NUM_CHANS; ++i)
      e.r[i] = targetPos[i] / 100.0;
    e.moveT = moveTime;
    e.postT = holdTime;

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
