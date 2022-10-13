#include "util/util.hpp"

using namespace std::chrono_literals;

int main(int argc, char *argv[]) {
  // good example on using chrono correctly:
  // https://stackoverflow.com/a/67805863 another example on clever casting
  // from durations to numerical types:
  // https://serveanswer.com/questions/how-to-cast-a-double-into-std-chrono-milliseconds

  Eigen::Vector3d x, z, y, r, u;
  r << 1, 2, 0;

  // iterator-based indexing of Eigen matrices
  // for (auto ri : r)
  //     std::cout << ri << "\n";
  // index-based indexing of Eigen matrices
  // for (std::size_t ind = 0; ind != r.rows(); ++ind) {
  //     std::cout << r(ind) << "\n";
  // }

  Eigen::MatrixXd P0 = Eigen::MatrixXd::Identity(5, 5);
  Eigen::Matrix<double, 5, 1> xhat0; // current state, state estimate vectors

  // array to eigen matrix
  // int16_t u_saturated[5] = {1,2,3,4,5};
  // xhat0 << 0,
  //          0,
  //          0,
  //          1,
  //          1;
  // Eigen::MatrixXd usat = Eigen::Map<Eigen::Matrix<int16_t,5,1>
  // >(u_saturated).cast<double>(); Eigen::MatrixXd usat_ =
  // P0*Eigen::Map<Eigen::Matrix<int16_t,5,1> >(u_saturated).cast<double>();
  // std::cout << usat << "\n";
  // std::cout << usat_ << "\n";
  // std::cout << P0 << "\n";
  // std::cout << xhat0 << "\n";

  Eigen::VectorXd a, b, c;
  a.resize(3);
  b.resize(3);
  c.resize(5);
  a << 1, 2, 3;
  b << 1, 2, 4;
  c << 1, 2, 3, 4, 5;

  c(0) = a(0);
  c(1) = a(1);
  c(2) = a(2);
  a = b;
  info("a: {}", a);
  info("c: {}", c);

  info("slice of c: {}", c(b.array(), Eigen::all));

  Eigen::MatrixXd dynamicMatrix; // = Eigen::Matrix<double, 3, 3>({1,2,3},{4,5,6},{7,8,9});
  dynamicMatrix.resize(3, 3);
  dynamicMatrix << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  info("dynamicMatrix:\n{}", dynamicMatrix);

  Eigen::Matrix<double, 3, 3> fixedMatrix = dynamicMatrix;
  info("fixedMatrix:\n{}", fixedMatrix);

  Eigen::Matrix<int, 3, 1> selectedChs = Eigen::Matrix<int, 3, 1>(1, 2, 3);
  info("selectedChs:\n{}", selectedChs);

  Eigen::VectorXi prbsOrder = Eigen::Vector3i(0, 1, 2);



  return 0;
}

// double Cd_data[3][3] = {{2.9187, -0.7636, -1.2614},
//                         {-0.2795, 0.6576, -0.6295},
//                         {-0.2366, 0.0339, 0.1116}};
// cv::Mat Cd(3, 3, CV_64FC1, Cd_data);
// double K1_data[3][3] = {{3.0931, -0.8994, -1.1502},
//                         {-0.4497, 0.8751, -0.5475},
//                         {-0.5751, -0.5475, 1.3220}};
// cv::Mat K1(3, 3, CV_64FC1, K1_data);
// double K2_data[3][3] = {{-0.8588, 0.1509, 0.2972},
//                         {0.1361, -0.6018, 0.7731},
//                         {0.2925, 0.7572, 0.5554}};
// cv::Mat K2(3, 3, CV_64FC1, K2_data);
// std::cout << "rows: " << Cd.rows << ", cols: " << Cd.cols << "\n";
// std::cout << Cd << "\n";
// std::cout << K1 << "\n";
// std::cout << K2 << "\n";

/*
        // double K1[3] = {-2.3510, 2.0327, -0.8836};    // state feedback gain
   (ch2)
        // double K1[3] = {-0.5290, -1.3344, 1.5552};    // state feedback gain
   (ch3) double K1[3] = {-3.1708, 0.6072, 0.2155};    // state feedback gain
   (ch1) double x = 0;                                 // state vector
        // double K2[3] = {1.3226, -1.1435, 0.4971};     // integral error
   feedback gain (ch2)
        // double K2[3] = {0.4218, 1.0641, -1.2403};     // integral error
   feedback gain (ch3) double K2[3] = {2.4159, -0.4626, -0.1642};     //
   integral error feedback gain (ch1) double z = 0; // integral error double y =
   0;                                 // measurement vector double r = 50; //
   reference input int16_t u[3] = {0,0,0};                       // control
   signal
*/
