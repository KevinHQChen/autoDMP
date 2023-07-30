#pragma once

#include "gui/components/button.hpp"
#include "gui/components/implot_helpers.hpp"
#include "gui/components/slider.hpp"
#include "window.hpp"

#include "ctrl/supervisor.hpp"

#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace gui {

namespace py = pybind11;
using namespace py::literals;

class CtrlWindow : public Window {
  Supervisor *sv_;
  ImProc *imProc_;

  int no;
  std::vector<double> rPixels;
  std::vector<int> yMax;

  // displaying/modifying events, states
  ImGuiTableFlags tableFlags = ImGuiTableFlags_BordersV | ImGuiTableFlags_BordersOuterH |
                               ImGuiTableFlags_Resizable | ImGuiTableFlags_RowBg |
                               ImGuiTableFlags_NoBordersInBody;
  const char *dragFmt = "%.01f", *txtFmt = "%.1f", *intFmt = "%d";

  event_bus currEv_, newEv_;
  int srcState, destState, moveTime, holdTime;
  double excitationAmp, paramLowerBound, covModification;
  double excitationAmpMin{0.0}, excitationAmpMax{10.0}, paramLowerBoundMin{0.0005},
      paramLowerBoundMax{0.1}, covModificationMin{0.0005}, covModificationMax{0.1};
  double uwt, umax;
  double uwtMin{0.0}, uwtMax{10.0}, umaxMin{80.0}, umaxMax{200.0};
  double minValMin{-10.0f}, minValMax{0.0f}, maxValMin{0.0f}, maxValMax{10.0f};

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

  // for sysID
  std::vector<double> timeVec_, uVec_;
  double minVal_, maxVal_;
  int order_;

  void renderAddEventDialog();
  void renderDropletGenDialog();
  void renderEventQueueContents();
  void renderSupervisorStatus();
  void renderControllerTuningDialog();
  void renderSysIdDialog();

  void loadEventsFromFile(const std::string &filename);

  void generateExcitationSignal(double minVal, double maxVal, int order);

public:
  bool ctrlSetupVisible_{false}, ctrlVisible_{false};

  CtrlWindow(Supervisor *sv, ImProc *imProc);
  ~CtrlWindow();
  void render() override;
};

} // namespace gui
