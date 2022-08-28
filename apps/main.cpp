#include <functional>
#include <optional>

// This file will be generated automatically when you run the CMake
// configuration step. It creates a namespace called `autoDMP`. You can modify
// the source template at `configured_files/config.hpp.in`. #include
// <internal_use_only/config.hpp>

#include "cam/cam.hpp"
#include "util/util.hpp"

using namespace std;
using namespace spdlog;

// NOLINTNEXTLINE(bugprone-exception-escape)
int main() {
  toml::value conf = toml::parse("config/setup.toml");
  info("Parsed setup.toml: {}", conf);
  auto &camConf = toml::find(conf, "cam");

  /* VIDEO SOURCE SETUP */
  cam *onlineCam;
  cv::VideoCapture *offlineCam;
  if (toml::get<std::string>(camConf["source"]) == "Andor") {
    onlineCam = new cam(0, conf);
    // start camera (allocate circular buffer to store frames)
    // timerInterval sets the size of buffers needed (min size of 1)
    // multiplied by the framerate(pretty sure)
    // if timerInterval is less than 1000ms, only 1 buffer (i.e .40 frames) is
    // needed
    onlineCam->start((int)(100 / 1000)); // timerInterval of 100ms
  } else if (toml::get<std::string>(camConf["source"]) == "File")
    offlineCam = new cv::VideoCapture(toml::get<std::string>(camConf["File"]));
  else if (toml::get<std::string>(camConf["source"]) == "Webcam")
    offlineCam = new cv::VideoCapture(0);

  // bool run = false;

  // // store command line arguments in config object
  // config conf;
  // conf.parse(argc, argv);

  // std::optional<std::string> message;
  // app.add_option("-m,--message", message, "A message to print back
  // out");

  // CLI11_PARSE(app, argc, argv);

  // if (show_version) {
  //   // fmt::print("{}\n", autoDMP::cmake::project_version);
  //   fmt::print("{}\n", "0.0.1");
  //   return EXIT_SUCCESS;
  // }

  // // Use the default logger (stdout, multi-threaded, colored)
  // spdlog::info("Hello, {}!", "World");

  // cv::namedWindow("Display window");
  // cv::Mat image;
  // char key;
  // cv::VideoCapture* cam = new cv::VideoCapture(0);

  // if (!cam->isOpened()) {
  //   fmt::print("cannot open camera\n");
  // }

  // while(true) {
  //   if(!cam->read(image)) {
  //     fmt::print("cannot read image\n");
  //     return 1;
  //   } else {
  //     cv::imshow("Display window", image);
  //     key = (char) cv::waitKey(1);
  //     if(key != -1) {
  //       cv::destroyAllWindows();
  //       break;
  //     }
  //   }
  // }

  // if (message) {
  //   fmt::print("Message: '{}'\n", *message);
  // } else {
  //   fmt::print("No Message Provided :( (use -m, --message then provide
  //   a message.)\n");
  // }
  // } catch (const toml::parse_error &e) {
  //   error("Unhandled exception in main: {}", e.what());
  // }
}
