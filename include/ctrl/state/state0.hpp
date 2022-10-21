#pragma once

#include "ctrl/state/state.hpp"

class State0 : public State {
public:
  // system matrices
  Vector1ui ch;
  Vector1d Ad, Ad_, Bd, Cd, Cd_, CdInv, K1, K2, Qw, Rv;
  // dynamic (changes based on observer estimation error)
  Vector1d P0, P, Ko, temp, tempInv;

  State0(Supervisor *supervisor);
  ~State0();

  virtual bool measurementAvailable() override;
  virtual void updateMeasurement() override;

  // update instantaneous trajectory vectors
  // (only called when new measurements are available)
  template <int dim> void updateMeasurement() {
    for (int i = 0; i != ch<dim>.rows(); ++i) {
      // update measurement vectors dy, y
      if (trueMeasAvail[ch<dim>(i)])
        dy(ch<dim>(i)) = sv_->imProc->procDataQArr[ch<dim>(i)]->get().y - yref(ch<dim>(i));
      else if (stateTransitionCondition) // provide constant output disturbance toward
                                         // junction
        dy(ch<dim>(i)) += dyref(ch<dim>(i));
      else
        dy(ch<dim>(i)) = dyhat(ch<dim>(i));
      y(ch<dim>(i)) = dy(ch<dim>(i)) + yref(ch<dim>(i));

      // update time step dt for numerical integration
      if (!firstMeasAvail[ch<dim>(i)])
        firstMeasAvail[ch<dim>(i)] = true;
      else // measure the time difference between consecutive measurements
        dt[ch<dim>(i)] = steady_clock::now() - prevCtrlTime[ch<dim>(i)];
      prevCtrlTime[ch<dim>(i)] = steady_clock::now();

      // update integral error based on time step for each channel's data
      z(ch<dim>(i)) += -dy(ch<dim>(i)) * dt[ch<dim>(i)].count();
    }
  }

  /*
   * @brief: performs the required tasks for state transition corresponding to the received event
   * @param: event
   */
  virtual void handleEvent(Event *event) override;
};
