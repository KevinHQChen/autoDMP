#pragma once

#include "ctrl/state/state.hpp"

class State1 : public State {
public:
  State1(Supervisor *supervisor, Eigen::Vector3d uref_);
  ~State1();

  // update instantaneous trajectory vectors
  // (only called when new measurements are available)
  template <int dim> void updateMeasurement() {
    for (int i = 0; i != ch<dim>.rows(); ++i) {
      // update measurement vectors dy, y
      if (trueMeasAvail[ch<dim>(i)])
        dy(ch<dim>(i)) = sv_->imProc->procDataQArr[ch<dim>(i)]->get().y - yref(ch<dim>(i));
      else if (stateTransitionCondition) // assume interface is stuck at junction (i.e. yref)
        dy(ch<dim>(i)) = yref0(ch<dim>(i)) - yref(ch<dim>(i));
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
