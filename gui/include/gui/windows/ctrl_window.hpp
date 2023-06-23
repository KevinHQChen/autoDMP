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
  const char *dragFmt = "%.01f", *txtFmt = "%.1f";

  event_bus currEv_, newEv_;
  int srcState, destState, moveTime, holdTime;
  float excitationAmp, paramLowerBound, covModification;

  int openAction = -1;
  double dropletLength = 0, dNeck = 0, dPlug = 0, wCh = 0;
  int numDroplets = 0;

  event_bus getEvent(std::vector<double> r, int preT, int moveT, int postT) {
    event_bus e;
    for (int i = 0; i < 2 * sv_->no; ++i)
      e.r[i] = r[i];
    e.preT = preT;
    e.moveT = moveT;
    e.postT = postT;

    return e;
  }

  void renderAddEventDialog();
  void renderDropletGenDialog();
  void renderEventQueueContents();
  void renderSupervisorStatus();
  void renderControllerTuningDialog();

public:
  bool ctrlSetupVisible_{false}, ctrlVisible_{false};

  CtrlWindow(Supervisor *sv);
  ~CtrlWindow();
  void render() override;
};

} // namespace gui
