#include "ctrl/supervisor.hpp"
#include "ctrl/state/state.hpp"
#include "ctrl/state/state0.hpp"
// #include "ctrl/state/state1.hpp"

State0::State0(Supervisor* sv) : State(sv) {
    yref = yref_traj[0];
}

State0::~State0() {
    // clean up any resources used by current state here
}

Eigen::Matrix<int16_t, 3, 1> State0::step() {
    ctrlSigs = getCtrlSigs(model, sv_->currPos, traj->currTraj);

    // stabilize at initial value of State0 trajectory
    if (abs(sv_->currPos - currTraj) > threshold)
        return ctrlSigs;




    if (currPos[0] )
        supervisor_->updateState<State1>();

    return ctrlSigs;
}
