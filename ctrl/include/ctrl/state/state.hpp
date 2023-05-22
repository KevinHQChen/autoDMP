#pragma once

#include "ctrl/supervisor.hpp"
#include "util/util.hpp"

// pump layout
//   top-down view
//     P3      P4
//     oil     outlet
//     | ch2   | ch3
//      \     /
//        \ /
//         |
//         | ch1
//         |
//        / \
//      /     \
//     |       |
//     P2      P1
//     water  water

//   camera view
//      y(0)
//      u(0)
//      water
//      P1 + P2
//         |
//         | ch1
//         |
//        / \
//      /     \
//     | ch2   | ch3
//     oil     outlet
//     P3      P4
//     u(1)    u(2)
//     y(1)    y(2)

template <int state, int numX, int numY, int numU> struct MDL {
  Supervisor *sv; // TODO remove this dependency somehow

  Eigen::Matrix<double, numX, 1> dxhat{Eigen::Matrix<double, numX, 1>::Zero()};
  Eigen::Matrix<double, numX, numX> Ad, Ad_;
  Eigen::Matrix<double, numX, numU> Bd;
  Eigen::Matrix<double, numY, numX> Cd;
  Eigen::Matrix<double, numX, numY> Cd_;
  Eigen::Matrix<double, numU, numX> K1;
  Eigen::Matrix<double, numY, numY> K2, Rv;
  Eigen::Matrix<double, numX, numX> Qw;
  // dynamic (changes based on observer estimation error)
  Eigen::Matrix<double, numX, numX> P0{Eigen::Matrix<double, numX, numX>::Identity()}, P{P0};
  Eigen::Matrix<double, numX, numY> Ko;
  Eigen::Matrix<double, numY, numY> temp, tempInv;

  MDL(Supervisor *sv)
      : Ad(openData(sv->getConfPath() + "state" + std::to_string(state) + "/Ad.txt")),
        Ad_(openData(sv->getConfPath() + "state" + std::to_string(state) + "/Ad_.txt")),
        Bd(openData(sv->getConfPath() + "state" + std::to_string(state) + "/Bd.txt")),
        Cd(openData(sv->getConfPath() + "state" + std::to_string(state) + "/Cd.txt")),
        Cd_(openData(sv->getConfPath() + "state" + std::to_string(state) + "/Cd_.txt")),
        K1(openData(sv->getConfPath() + "state" + std::to_string(state) + "/K1.txt")),
        K2(openData(sv->getConfPath() + "state" + std::to_string(state) + "/K2.txt")),
        Rv(openData(sv->getConfPath() + "state" + std::to_string(state) + "/Rv.txt")),
        Qw(openData(sv->getConfPath() + "state" + std::to_string(state) + "/Qw.txt")) {}
};

class State {
  template <bool COND, int A, int B> struct IF {
    enum { val = A };
  };

  template <int A, int B> struct IF<false, A, B> {
    enum { val = B };
  };

  // stores min of A and B
  template <int A, int B> struct MIN : IF < A<B, A, B> {};

protected:
  Supervisor *sv_;

public:
  // from https://stackoverflow.com/a/15484513
  template <typename T, int dim, int newDim>
  Eigen::Matrix<T, newDim, 1> to(Eigen::Matrix<T, dim, 1> p) {
    Eigen::Matrix<int, newDim, 1> newp = Eigen::Matrix<T, newDim, 1>::Zero();

    newp.template head<MIN<dim, newDim>::val>() = p.template head<MIN<dim, newDim>::val>();

    return newp;
  }

  // instantaneous trajectory vectors
  // control signal vectors
  Eigen::Vector3d du{Eigen::Vector3d::Zero()}, uref, u, usat;
  // integral error vectors: z = int (r - y) dt
  Eigen::Vector3d z0{Eigen::Vector3d::Zero()}, z{z0};
  time_point<steady_clock> initTime{steady_clock::now()};
  time_point<steady_clock> prevCtrlTime[3] = {initTime, initTime, initTime};
  duration<double> dt[3] = {0s, 0s, 0s}; // this somehow defaults to milliseconds
  // state/output vectors
  Eigen::Vector3d yrefScale;
  Eigen::Vector3d y{Eigen::Vector3d::Zero()}, yref0, yref, ytrans{Eigen::Vector3d::Zero()};
  Eigen::Vector3d dy, dyref{Eigen::Vector3d::Zero()};
  Eigen::Vector3d yhat, dyhat{Eigen::Vector3d::Zero()}, ytilde;
  Eigen::Vector3d yDest;

  bool firstMeasAvail[3] = {false, false, false};
  bool measAvail[3] = {true, true, true};
  bool trueMeasAvail[3] = {true, true, true};
  bool obsv[3] = {false, false, false};
  bool stateTransitionCondition = false;
  bool startEvent = false;

  State(Supervisor *sv, Eigen::Vector3d uref_, Eigen::Vector3d yrefScale);
  virtual ~State();

  // provide common interface for member functions using pure virtual methods
  // (this makes State an abstract class)

  // check for new measurements on selected channels
  virtual bool measurementAvailable() = 0;

  /* @brief: update instantaneous trajectory vectors (only called when new measurements are
   * available)
   */
  virtual void updateMeasurement() = 0;

  /*
   * @brief: performs the required tasks for state transition corresponding to the received event
   * @param: event
   */
  virtual void handleEvent(Event *event) = 0;

  /*
   * @brief: generates control signals at each time step (only called when new measurements are
   * available)
   */
  virtual Eigen::Vector3d step() = 0;

  template <int dim> bool measurementAvailable(Eigen::Matrix<unsigned int, dim, 1> &ch) {
    bool tmpMeasAvail = true;
    for (int i = 0; i != ch.rows(); ++i) {
      trueMeasAvail[ch(i)] = !sv_->p[ch(i)].found;
      measAvail[ch(i)] =
          duration_cast<milliseconds>(steady_clock::now() - prevCtrlTime[ch(i)]).count() >= 25;

      // only allow true measurements if we're within channel boundaries (35%->85%)
      // otherwise allow estimated measurements from the observer
      // TODO WE SHOULD PROBABLY USE Y INSTEAD OF YREF, LOOK INTO IT MORE
      tmpMeasAvail &=
          (yref(ch(i)) < 0.85 * yrefScale(ch(i)) && yref(ch(i)) > 0.35 * yrefScale(ch(i)))
              ? trueMeasAvail[ch(i)]
              : measAvail[ch(i)];
    }

    return tmpMeasAvail;
  }

  template <int dim> void updateMeasurement(Eigen::Matrix<unsigned int, dim, 1> &ch, int rot) {
    Pose p;
    for (int i = 0; i != ch.rows(); ++i) {
      // update measurement vectors dy, y
      if (trueMeasAvail[ch(i)])
        p = sv_->p[ch(i)];

      if (trueMeasAvail[ch(i)] && (rot == -1))
        dy(ch(i)) = p.loc[0].y - yref(ch(i));
      else if (stateTransitionCondition) { // assume interface is stuck at junction
        if (ytrans(ch(i)) < yref(ch(i)))
          ytrans(ch(i)) = yref(ch(i));
        dy(ch(i)) = ytrans(ch(i)) - yref(ch(i));
      } else // use estimated value from kalman observer
        dy(ch(i)) = dyhat(ch(i));
      y(ch(i)) = dy(ch(i)) + yref(ch(i));

      // TODO move everything after this point to the step function
      // update time step dt for numerical integration
      if (!firstMeasAvail[ch(i)])
        firstMeasAvail[ch(i)] = true;
      else // measure the time difference between consecutive measurements
        dt[ch(i)] = steady_clock::now() - prevCtrlTime[ch(i)];
      prevCtrlTime[ch(i)] = steady_clock::now();

      // update integral error based on time step for each channel's data
      z(ch(i)) += -dy(ch(i)) * dt[ch(i)].count();

      // update uref to minimize du
      uref(ch(i)) += 0.1 * du(ch(i)) * dt[ch(i)].count();
    }
  }

  template <int state, int numX, int numY, int numU>
  Eigen::Vector3d step(Eigen::Matrix<unsigned int, numY, 1> &ch,
                                    MDL<state, numX, numY, numU> *mdl) {
    // update state x and control signal u based on new measurements
    // Kalman observer
    // prediction
    // dxhat = dxhat_prev - dxref, Cd*dxref = dyref
    mdl->dxhat = mdl->Ad * mdl->dxhat + mdl->Bd * du(ch.array(), Eigen::all);
    mdl->P = mdl->Ad * mdl->P * mdl->Ad_ + mdl->Qw;
    // correction
    mdl->temp = mdl->Cd * mdl->P * mdl->Cd_ + mdl->Rv;
    mdl->tempInv = mdl->temp.inverse();
    mdl->Ko = mdl->P * mdl->Cd_ * mdl->tempInv;
    dyhat(ch.array(), Eigen::all) = mdl->Cd * mdl->dxhat;
    mdl->dxhat =
        mdl->dxhat + mdl->Ko * (dy(ch.array(), Eigen::all) - dyhat(ch.array(), Eigen::all));

    // apply updated control signals (with saturation limits) to pump
    // LQI control law
    du(ch.array(), Eigen::all) = -mdl->K1 * mdl->dxhat - mdl->K2 * z(ch.array(), Eigen::all);
    u = uref + du;

    // apply saturation (+/- 20)
    int sat = 40;
    for (int i = 0; i != du.rows(); ++i) {
      if (du(i) < -sat)
        du(i) = -sat;
      else if (du(i) > sat)
        du(i) = sat;
    }
    usat = uref + du;

    return usat;
  }
};
