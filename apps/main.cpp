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

int
main(int, const char**)
{
#if __cplusplus__ >= 202000ULL
    ImGuiWrapConfig config{.windowTitle_ = "dear::Example", .width_ = 1280, .height_ = 600};
#else
    ImGuiWrapConfig config{};
    config.windowTitle_ = "dear::Example";
    config.width_       = 800;
    config.height_      = 600;
#endif

    return imgui_main(config, windowFn);
}

#if 0
// NOLINTNEXTLINE(bugprone-exception-escape)
int main() {
#if 0
  GUIRenderer renderer;
  renderer.InitGUI();
  renderer.Render();
#endif

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

  GUIRenderer renderer;
  renderer.InitGUI();
  renderer.Render();

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
