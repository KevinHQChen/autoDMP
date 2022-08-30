#include <functional>
#include <optional>

// This file will be generated automatically when you run the CMake
// configuration step. It creates a namespace called `autoDMP`. You can modify
// the source template at `configured_files/config.hpp.in`. #include
// <internal_use_only/config.hpp>

#include "cam/cam.hpp"
#include "gui/gui.hpp"
#include "util/util.hpp"

using namespace std;
using namespace spdlog;

// NOLINTNEXTLINE(bugprone-exception-escape)
int main() {
  info("setup DSL: status code = {}",
       SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER));

  // glsl_version corresponds to OpenGL version (use glxinfo | grep version to
  // get OpenGL version) (see
  // https://www.khronos.org/opengl/wiki/Core_Language_(GLSL)#Version)
  const char *glsl_version = "#version 460";
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);

  // Create window with graphics context
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
  SDL_WindowFlags window_flags =
      (SDL_WindowFlags)(SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE |
                        SDL_WINDOW_ALLOW_HIGHDPI);
  SDL_Window *window = SDL_CreateWindow(
      "Dear ImGui SDL2+OpenGL3 example", SDL_WINDOWPOS_CENTERED,
      SDL_WINDOWPOS_CENTERED, 1280, 720, window_flags);
  SDL_GLContext gl_context = SDL_GL_CreateContext(window);
  SDL_GL_MakeCurrent(window, gl_context);
  SDL_GL_SetSwapInterval(1); // Enable vsync

  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  (void)io;
  io.ConfigFlags |=
      ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
  // io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad
  // Controls

  // Setup Dear ImGui style
  ImGui::StyleColorsDark();
  // ImGui::StyleColorsLight();

  // Setup Platform/Renderer backends
  ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
  ImGui_ImplOpenGL3_Init(glsl_version);

  // Our state
  bool show_demo_window = true;
  bool show_another_window = false;
  ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

  bool done = false;
  while (!done) {
    // Poll and handle events (inputs, window resize, etc.)
    // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to
    // tell if dear imgui wants to use your inputs.
    // - When io.WantCaptureMouse is true, do not dispatch mouse input data to
    // your main application, or clear/overwrite your copy of the mouse data.
    // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input
    // data to your main application, or clear/overwrite your copy of the
    // keyboard data. Generally you may always pass all inputs to dear imgui,
    // and hide them from your application based on those two flags.
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      ImGui_ImplSDL2_ProcessEvent(&event);
      if (event.type == SDL_QUIT)
        done = true;
      if (event.type == SDL_WINDOWEVENT &&
          event.window.event == SDL_WINDOWEVENT_CLOSE &&
          event.window.windowID == SDL_GetWindowID(window))
        done = true;
    }

    // Start the Dear ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplSDL2_NewFrame();
    ImGui::NewFrame();

    // 1. Show the big demo window (Most of the sample code is in
    // ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear
    // ImGui!).
    if (show_demo_window)
      ImGui::ShowDemoWindow(&show_demo_window);

    // 2. Show a simple window that we create ourselves. We use a Begin/End pair
    // to created a named window.
    {
      static float f = 0.0f;
      static int counter = 0;

      ImGui::Begin("Hello, world!"); // Create a window called "Hello, world!"
                                     // and append into it.

      ImGui::Text("This is some useful text."); // Display some text (you can
                                                // use a format strings too)
      ImGui::Checkbox(
          "Demo Window",
          &show_demo_window); // Edit bools storing our window open/close state
      ImGui::Checkbox("Another Window", &show_another_window);

      ImGui::SliderFloat("float", &f, 0.0f,
                         1.0f); // Edit 1 float using a slider from 0.0f to 1.0f
      ImGui::ColorEdit3(
          "clear color",
          (float *)&clear_color); // Edit 3 floats representing a color

      if (ImGui::Button("Button")) // Buttons return true when clicked (most
                                   // widgets return true when edited/activated)
        counter++;
      ImGui::SameLine();
      ImGui::Text("counter = %d", counter);

      ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                  1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
      ImGui::End();
    }

    // 3. Show another simple window.
    if (show_another_window) {
      ImGui::Begin(
          "Another Window",
          &show_another_window); // Pass a pointer to our bool variable (the
                                 // window will have a closing button that will
                                 // clear the bool when clicked)
      ImGui::Text("Hello from another window!");
      if (ImGui::Button("Close Me"))
        show_another_window = false;
      ImGui::End();
    }

    // Rendering
    ImGui::Render();
    glViewport(0, 0, (int)io.DisplaySize.x, (int)io.DisplaySize.y);
    glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w,
                 clear_color.z * clear_color.w, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    SDL_GL_SwapWindow(window);
  }

  // Cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplSDL2_Shutdown();
  ImGui::DestroyContext();

  SDL_GL_DeleteContext(gl_context);
  SDL_DestroyWindow(window);
  SDL_Quit();

  /*
  ordered_value conf = toml::parse<toml::discard_comments, tsl::ordered_map>(
      "config/setup.toml");
  info("Config type: {}", type_name<decltype(conf)>());
  info("Parsed config: {}", conf);

  auto &camConf = toml::find(conf, "cam");
  info("Camconf type: {}", type_name<decltype(camConf)>());
  info("Parsed camConf: {}", camConf);

  // VIDEO SOURCE SETUP
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

  info("Starting image capture...");
  cv::Mat currentImage = cv::Mat(0, 0, CV_16UC1);
  cv::Mat currentImageRaw = cv::Mat(0, 0, CV_16UC1);
  cv::namedWindow("Display window");
  bool imCapSuccess;
  char key;

  while (true) {
    imCapSuccess = !onlineCam->process(currentImage);
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

  delete onlineCam;
  */

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
