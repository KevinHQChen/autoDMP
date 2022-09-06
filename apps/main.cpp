// config.hpp will be generated automatically when you run the CMake
// configuration step. It creates a namespace called `autoDMP`. You can modify
// the source template at `configured_files/config.hpp.in`. #include
// <internal_use_only/config.hpp>

#include "gui/gui.hpp"

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int, const char **) {
  GUI *gui = new GUI();
  gui->startGUIThread();
  // ordered_value conf = toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml");
  // Cam *cam = new Cam(0, conf);
  // cv::Mat currentImage = cv::Mat(0, 0, CV_16UC1);
  // cv::namedWindow("Display window");
  // bool imCapSuccess;
  // char key;

  // cam->start((int)(100 / 1000)); // timerInterval of 100ms
  // while (true) {
  //   imCapSuccess = cam->process(currentImage);
  //   if (currentImage.type() == CV_16UC1)
  //     info("16 bit image");
  //   if (!currentImage.empty())
  //     info("Image not empty");

  //   if (imCapSuccess) {
  //     cv::imshow("Display window", currentImage);
  //     key = (char)cv::waitKey(1);
  //     if (key != -1) {
  //       cv::destroyAllWindows();
  //       break;
  //     }
  //   } else {
  //     error("cannot read image");
  //   }
  // }
  // delete cam;
}
