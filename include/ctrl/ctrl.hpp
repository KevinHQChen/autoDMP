#pragma once

#include "util/util.hpp"
#include "improc/improc.hpp"
#include "pump/pump.hpp"
#include "ctrl/supervisor.hpp"


class adgModel {
        public:
                // ADG 0: move ch0 interface through junction until it appears in ch1/ch2
                //     selCh = 0
                //     uref = [60, 40, 60]
                //     yref_traj:
                //         ch0: ramp (0.5->1.05) over 30s, scale by channel height
                // ADG 1: move ch2 interface to desired length while holding ch1 in place
                //     selCh = [1, 2]
                //     uref = [65, 45, 65]
                //     yref_traj:
                //         ch1: constant (0.8), scale by channel height
                //         ch2: constant (0.9), scale by channel height
                // ADG 2: move ch1 interface toward junction while holding ch2 until split occurs
                //     selCh = [1, 2]
                //     uref = [65, 45, 65]
                //     yref_traj:
                //         ch1: ramp (0.8->1.05) over 15s, scale by channel height
                //         ch2: constant (0.9), scale by channel height
                // ADG 3: FILL IN
                //     selCh = [0, 2]
                //     uref = [?, ?, ?]
                //     yref_traj:
                //         ch0: FILL IN
                //         ch2: FILL IN
                adgModel(int ADG, ChannelPose chanPose);
                ~adgModel();

                // vector of selected channels
                Eigen::VectorXi selCh;
                int numCh, freq;

                // system matrices (todo: make matrices fixed-size for constants)
                // static for each ADG step
                Eigen::MatrixXd Ad, Ad_, Bd, Cd, Cd_, CdInv, K1, K2, Qw, Rv;
                // dynamic (changes based on observer estimation error)
                Eigen::MatrixXd P0, Ko, temp, tempInv;

                // instantaneous trajectory vectors
                //     current integral error, measurement vectors
                //     resize for different model structures, e.g. z.resize(3);
                Eigen::VectorXd z0; // integral error: z = int (r - y) dt
                Eigen::VectorXd y, dy, dyhat, dxhat;
                Eigen::VectorXd yref0, dyref, yref_scale;
                Eigen::Vector3d uref; // control signal vectors
                // obsv error: ytilde = (y - yhat)
                // (by default all vectors are column vectors, so we good)

                // reference trajectory
                // ADG 0:
                //     0->1.1 ramp over 30s, scale by channel height/2, offset by channel height/2
                // ADG 1:
                //
                Eigen::MatrixXd yref_traj;
        private:
                int _ADG;
};

void ctrl(config conf,
          std::vector<QueueFPS<cv::Point>*>& procDataQueues,
          std::vector<QueueFPS<cv::Mat>*>& procFramesQueues,
          ChannelPose& chanPose,
          QueueFPS<uint16_t>& ctrlDataQueue, // placeholder for shared memory idea for interfacing to nodejs
          bool& run);

void sendCtrlSignals(config conf, pump* pp, Eigen::Matrix<int16_t, 3, 1> usat);

void saveCtrlData(QueueFPS<uint16_t>& ctrlDataQueue, Eigen::VectorXd vect);
