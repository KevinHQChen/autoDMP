#pragma once

#include "ctrl/state/state.h"

class State0 : State {
    static Model model;
    static Eigen::MatrixXd yref_traj;
    double yref;

    int rotAngle;
    cv::Rect chanBBox;

    public:
        State0(Supervisor* supervisor);
        ~State0();
        virtual Eigen::Matrix<int16_t, 3, 1> step() override;
};
