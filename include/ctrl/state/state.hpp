#pragma once

#include "ctrl/supervisor.hpp"
#include "util/util.hpp"

// pump layout
//   top-down view
//     P3      P4
//     oil     water
//     | ch2   | ch1
//      \     /
//        \ /
//         |
//         | ch3
//         |
//        / \
//      /     \
//     |       |
//     P2      P1
//     outlet  outlet

//   camera view
//      y(2)
//      u(2)
//      outlet
//      P1 + P2
//         |
//         | ch3
//         |
//        / \
//      /     \
//     | ch2   | ch1
//     oil     water
//     P3      P4
//     u(1)    u(0)
//     y(1)    y(0)

struct Event;
class Supervisor; // forward declaration

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
  // system matrices
  template <int dim>
  static Eigen::Matrix<double, dim, dim>
      // constant (for each state)
      Ad, Ad_, Bd, Cd, Cd_, CdInv, K1, K2, Qw, Rv,
      // variable (changes based on observer estimation error)
      P0, P, Ko, temp, tempInv;

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
  Eigen::Vector3d y{Eigen::Vector3d::Zero()}, yref0, yref;
  Eigen::Vector3d dy, dyref;
  Eigen::Vector3d yhat, dxhat{Eigen::Vector3d::Zero()}, dyhat{Eigen::Vector3d::Zero()}, ytilde;
  Eigen::Vector3d yDest;

  bool firstMeasAvail[3] = {false, false, false};
  bool measAvail[3] = {true, true, true};
  bool trueMeasAvail[3] = {true, true, true};
  bool obsv[3] = {false, false, false};
  bool stateTransitionCondition = false;
  bool startEvent = false;

  // sysID specific
  int stp = 0;

  // from https://stackoverflow.com/a/15484513
  template <typename T, int dim, int newDim>
  Eigen::Matrix<T, newDim, 1> to(Eigen::Matrix<T, dim, 1> p) {
    Eigen::Matrix<int, newDim, 1> newp = Eigen::Matrix<T, newDim, 1>::Zero();

    newp.template head<MIN<dim, newDim>::val>() = p.template head<MIN<dim, newDim>::val>();

    return newp;
  }

  State(Supervisor *sv, Eigen::Vector3d uref_, Eigen::Vector3d yrefScale);
  virtual ~State();

  // check for new measurements on selected channels
  template <int dim> bool measurementAvailable(Eigen::Matrix<unsigned int, dim, 1> ch) {
    bool tmpMeasAvail = true;
    for (int i = 0; i != ch.rows(); ++i) {
      trueMeasAvail[ch(i)] = !sv_->imProc->procDataQArr[ch(i)]->empty();
      measAvail[ch(i)] =
          duration_cast<milliseconds>(steady_clock::now() - prevCtrlTime[ch(i)]).count() >= 25;

      if (obsv[ch(i)] || (!obsv[ch(i)] && trueMeasAvail[ch(i)] && !stateTransitionCondition)) {
        tmpMeasAvail &= trueMeasAvail[ch(i)];
        obsv[ch(i)] = true;
      } else
        tmpMeasAvail &= measAvail[ch(i)];
    }

    return tmpMeasAvail;
  }

  virtual bool measurementAvailable() = 0;
  virtual void updateMeasurement() = 0;

  /*
   * @brief: generates control signals at each time step (only called when new measurements are
   * available)
   */
  template <int dim> Eigen::Matrix<int16_t, 3, 1> step() {
    // update state x and control signal u based on new measurements
    // Kalman observer
    // prediction
    // dxhat = dxhat_prev - dxref, Cd*dxref = dyref
    dxhat(ch<dim>.array(), Eigen::all) = Ad<dim> * dxhat(ch<dim>.array(), Eigen::all) +
                                         Bd<dim> * du(ch<dim>.array(), Eigen::all) -
                                         CdInv<dim> * dyref(ch<dim>.array(), Eigen::all);
    P<dim> = Ad<dim> * P<dim> * Ad_<dim> + Qw<dim>;
    // correction
    temp<dim> = Cd<dim> * P<dim> * Cd_<dim> + Rv<dim>;
    tempInv<dim> = temp<dim>.inverse();
    Ko<dim> = P<dim> * Cd_<dim> * tempInv<dim>;
    dyhat(ch<dim>.array(), Eigen::all) = Cd<dim> * dxhat(ch<dim>.array(), Eigen::all);
    dxhat(ch<dim>.array(), Eigen::all) =
        dxhat(ch<dim>.array(), Eigen::all) +
        Ko<dim> * (dy(ch<dim>.array(), Eigen::all) - dyhat(ch<dim>.array(), Eigen::all));

    // apply updated control signals (with saturation limits) to pump
    // LQI control law
    du(ch<dim>.array(), Eigen::all) =
        -K1<dim> * dxhat(ch<dim>.array(), Eigen::all) - K2<dim> * z(ch<dim>.array(), Eigen::all);
    u = uref + du;

    // apply saturation (+/- 20)
    for (int i = 0; i != du.rows(); ++i) {
      if (du(i) < -40)
        du(i) = -40;
      else if (du(i) > 40)
        du(i) = 40;
    }
    usat = uref + du;

    return usat.cast<int16_t>();
  }

  virtual void handleEvent(Event *event) = 0;
};
