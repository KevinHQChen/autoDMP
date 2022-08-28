#include <functional>
#include <iostream>
#include <optional>

#include <CLI/CLI.hpp>
#include <spdlog/spdlog.h>

// This file will be generated automatically when you run the CMake
// configuration step. It creates a namespace called `autoDMP`. You can modify
// the source template at `configured_files/config.hpp.in`. #include
// <internal_use_only/config.hpp>

#include "util/util.hpp"
#include <opencv2/opencv.hpp>

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, const char **argv) {
  try {
    CLI::App app;

    int item;

    app.add_flag("--simple");
    app.add_option("--item", item);
    app.set_config("--config", "config.toml", "Read config toml file", true);

    CLI11_PARSE(app, argc, argv);

    std::cout << app.config_to_str(true, true) << std::endl;

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
  } catch (const std::exception &e) {
    spdlog::error("Unhandled exception in main: {}", e.what());
  }
}
