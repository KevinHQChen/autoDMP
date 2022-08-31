// config.hpp will be generated automatically when you run the CMake
// configuration step. It creates a namespace called `autoDMP`. You can modify
// the source template at `configured_files/config.hpp.in`. #include
// <internal_use_only/config.hpp>

#include "cam/cam.hpp"
#include "gui/gui.hpp"
#include "util/util.hpp"

ImGuiWrapperReturnType
windowFn()
{
    static bool show_visualizer{true};
    static bool selected{false};

    // Returning a value will translate to an exit code.
    if (!show_visualizer)
        return 0;

    static bool is_quitting;
    dear::Begin("Visualizer", &show_visualizer) && []() {
        dear::MainMenuBar() && []() {
            dear::Menu("File") && []() {
                ImGui::MenuItem("Wibble", nullptr, &selected);
                ImGui::MenuItem("BOOM");
            };
            dear::Menu("Quitters") && []() {
                is_quitting = ImGui::MenuItem("QUIT NAOW");
            };
        };
        ImGui::Text("hello");
        dear::TabBar("##TabBar") && []() {
            dear::TabItem("Files") && []() {
                ImGui::Text("...files...");
            };
            dear::TabItem("Blueprints") && []() {
                ImGui::Text("...blueprints...");
            };
            dear::TabItem("Enums") && []() {
                ImGui::Text("...enums...");
            };
            dear::TabItem("Prototypes") && []() {
                ImGui::Text("...prototypes...");
            };
            dear::TabItem("More") && []() {
                dear::Child("hello") && []() {
                    dear::Group() && []() {
                        dear::Combo("me", "you") && []() {
                            dear::Tooltip() && []() {
                                ImGui::SetTooltip("You are now viewing a tooltip...");
                            };
                        };
                        dear::ListBox("list of things");
                    };
                };
            };
        };
    };
    if (is_quitting)
        return 0;

    // Returns "no value" (see std::optional)
    return {};
}

ImGuiWrapperReturnType
liveWebcam()
{
  static bool videoSourceSetup{true};
  static cv::Mat currentImage, dispImage;
  static Cam *cam;
  static bool updateTexture{false};
  static ordered_value conf = toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml");

  static bool selected{false};
  static bool is_quitting{false};

  if (videoSourceSetup) {
    cam = new Cam(0, conf);
    // start camera (allocate circular buffer to store frames)
    // timerInterval sets the size of buffers needed (min size of 1)
    // multiplied by the framerate(pretty sure)
    // if timerInterval is less than 1000ms, only 1 buffer (i.e .40 frames) is
    // needed
    cam->start((int)(100 / 1000)); // timerInterval of 100ms
    info("Starting image capture...");
    // currentImage = cv::Mat(0, 0, CV_16UC1);
    videoSourceSetup = false;
  }

  if (cam->process(currentImage)) {
    info("New image captured...");
    info("Image size: {} x {}", currentImage.rows, currentImage.cols);
    updateTexture = true;
    if (is_quitting) {
      cv::destroyAllWindows();
      delete cam;
      return 0;
    }
  } else {
    error("cannot read image\n");
    updateTexture= false;
  }

  dear::Begin("Webcam") && []() {
    ImGui::Text("Description");
    dear::MainMenuBar() && []() {
      dear::Menu("File") && []() {
        ImGui::MenuItem("Test", nullptr, &selected);
        is_quitting = ImGui::MenuItem("Quit");
      };
    };

    if (updateTexture && !currentImage.empty()) {
      GLuint textureID;
      // std::string filename = "./bird.jpeg";
      // cv::Mat resized_image;
      // cv::Mat image = cv::imread(filename);
      // int resized_width = 640;
      // double scale = static_cast<float>(resized_width) / currentImage.size().width;
      // info("scale: {}", scale);
      // cv::resize(currentImage, resized_image, cv::Size(0, 0), scale, scale);
      info("Updating texture...");
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
      // if (textureID == -1)
      glGenTextures(1, &textureID);
      glBindTexture(GL_TEXTURE_2D, textureID);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      // glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
      glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
      // glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
      // glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
      // glTexImage2D(GL_TEXTURE_2D, 0, GL_SINGLE_COLOR, resized_image.cols, resized_image.rows, 0, GL_SINGLE_COLOR,
      //              GL_UNSIGNED_BYTE, resized_image.ptr());
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, currentImage.cols, currentImage.rows, 0, GL_BGR, GL_UNSIGNED_BYTE, currentImage.ptr());
      glGenerateMipmap(GL_TEXTURE_2D);
      ImGui::Image((void *)(intptr_t)textureID, ImVec2(currentImage.cols, currentImage.rows));
      updateTexture = false;
    }
  };

  return {};
}

#if 1
int
main(int, const char**)
{
  ordered_value conf = toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml");
  info("Config type: {}", type_name<decltype(conf)>());
  info("Parsed config: {}", conf);

  // // VIDEO SOURCE SETUP
  // Cam *cam = new Cam(0, conf);
  // cam->start((int)(100 / 1000));
  // info("Starting image capture...");

  // GUIRenderer renderer;
  // renderer.InitGUI();
  // renderer.Render(cam);


  return imgui_main(toml::find(conf, "gui"), liveWebcam);
}
#endif

#if 0
// NOLINTNEXTLINE(bugprone-exception-escape)
int main() {
#if 0
  ordered_value conf = toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml");
  info("Config type: {}", type_name<decltype(conf)>());
  info("Parsed config: {}", conf);

  // VIDEO SOURCE SETUP
  Cam *cam = new Cam(0, conf);
  // start camera (allocate circular buffer to store frames)
  // timerInterval sets the size of buffers needed (min size of 1)
  // multiplied by the framerate(pretty sure)
  // if timerInterval is less than 1000ms, only 1 buffer (i.e .40 frames) is
  // needed
  cam->start((int)(100 / 1000)); // timerInterval of 100ms
  info("Starting image capture...");
  cv::Mat currentImage = cv::Mat(0, 0, CV_16UC1);
  cv::namedWindow("Display window");
  bool imCapSuccess;
  char key;

  // GUI SETUP (replaces the ugliness below)
  // std::thread captureThread(imgui_main, toml::find(config, "gui"), windowFn);

  while (true) {
    imCapSuccess = cam->process(currentImage);
    if (imCapSuccess) {
      cv::imshow("Display window", currentImage);
      key = (char)cv::waitKey(1);
      if (key != -1) {
        cv::destroyAllWindows();
        break;
      }
    } else {
      error("cannot read image\n");
    }
  }
  delete cam;
#endif
}
#endif
