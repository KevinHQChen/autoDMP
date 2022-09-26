#include "gui/gui.hpp"
#include <cstdio>

static void glfw_error_callback(int error, const char *description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

// imgui_main initializes an ImGui openGL/glfw backend and then runs
// the passed std::function<std::optional<int>()> is called repeatedly until the std::optional it
// returns has a value, which is then returned as the exit code.
int GUI::imguiMain() {
  // Setup window
  glfwSetErrorCallback(glfw_error_callback);
  if (glfwInit() == 0) {
    return 1;
  }

  // Decide GL+GLSL versions
#if defined(IMGUI_IMPL_OPENGL_ES2)
  // GL ES 2.0 + GLSL 100
  const char *glsl_version = "#version 100";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#else
  // GL 4.6 + GLSL 460
  // glsl_version corresponds to OpenGL version
  // (use glxinfo | grep version to get OpenGL version)
  // (see https://www.khronos.org/opengl/wiki/Core_Language_(GLSL)#Version)
  const char *glsl_version = "#version 460";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // 3.2+ only
#endif

  // Create window with graphics context
  window = glfwCreateWindow(guiConf.width, guiConf.height, guiConf.windowTitle.c_str(), nullptr,
                            nullptr);
  if (window == nullptr) {
    return 1;
  }

  glfwMakeContextCurrent(window);
  glfwSwapInterval(guiConf.enableVsync); // Enable vsync

  // Initialize OpenGL loader
  bool err = glewInit() != GLEW_OK;
  if (err) {
    fprintf(stderr, "Failed to initialize OpenGL loader!\n");
    return 1;
  }

  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  (void)io;
  if (guiConf.keyboardNav)
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
  io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;       // Enable Docking

  guiConf.startDark ? ImGui::StyleColorsDark() : ImGui::StyleColorsLight();

  // Setup Platform/Renderer backends
  /// TODO: Needs to be based on cmake config.
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(glsl_version);

  // io.Fonts->AddFontDefault();
  io.Fonts->AddFontFromFileTTF(guiConf.fontPath.c_str(), guiConf.fontSize);
  // ImGuiStyle& style = ImGui::GetStyle();
  // style.ScaleAllSizes(guiConf.scale);

  // Main loop
  const auto &clearColor = guiConf.clearColor;
  std::optional<int> exitCode{};

  while (!exitCode.has_value() && glfwWindowShouldClose(window) == 0) {
    // Poll and handle events (inputs, window resize, etc.)
    // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui
    // wants to use your inputs.
    // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main
    // application.
    // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main
    // application. Generally you may always pass all inputs to dear imgui, and hide them from
    // your application based on those two flags.
    glfwPollEvents();

    // Start the Dear ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    exitCode = this->render();

    // Rendering
    ImGui::Render();

    // NOLINTNEXTLINE(readability-isolate-declaration) input parameters to next call.
    int display_w{0}, display_h{0};
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);

    // setup the 'clear' background.
    glClearColor(clearColor[0] * clearColor[3], clearColor[1] * clearColor[3],
                 clearColor[2] * clearColor[3], clearColor[3]);
    glClear(GL_COLOR_BUFFER_BIT);

    // finialize the imgui render into draw data, and render it.
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    // swap the render/draw buffers so the user can see this frame.
    glfwSwapBuffers(window);

    // change the native (host) window size if requested.
    if (newSize.has_value()) {
      glfwSetWindowSize(window, newSize.value().first, newSize.value().second);
      newSize.reset();
    }
  }

  // Cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  glfwDestroyWindow(window);
  glfwTerminate();

  return exitCode.value_or(0);
}

void GUIFrame::updateTexture() {
  // mimic opencv grayscale image in opengl by making each color channel the same
  if (mat.empty()) {
    error("Image is empty");
    return;
  }

  // get image data type
  GLenum imgDataType;
  if (mat.type() == CV_8UC1) {
    imgDataType = GL_UNSIGNED_BYTE;
    // info("Image is CV_8UC1");
  } else if (mat.type() == CV_16UC1) {
    imgDataType = GL_UNSIGNED_SHORT;
    // info("Image is CV_16UC1");
  } else {
    error("Unsupported image data type");
    return;
  }

  cv::Mat tmp;
  cv::merge(std::vector<cv::Mat>{mat, mat, mat}, tmp);

  // https://github.com/ocornut/imgui/issues/4628
  if (texture == 0) {
    // update texture
    // glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    // create opengl texture identifier
    glGenTextures(1, &texture);
  }
  glBindTexture(GL_TEXTURE_2D, texture);
  // setup filtering parameters for display
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  // upload pixels into texture
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glPixelStorei(
      GL_UNPACK_ALIGNMENT,
      1); // https://stackoverflow.com/a/45486871
          // (https://www.khronos.org/opengl/wiki/Common_Mistakes#Texture_upload_and_pixel_reads)
  glTexImage2D(GL_TEXTURE_2D, // texture type
               0,             // pyramid level (for mip-mapping), 0 is the top level
               GL_RGB,        // internal color format to convert to
               // (https://www.khronos.org/opengl/wiki/Image_Format)
               tmp.cols, tmp.rows, // image width, height
               0,                  // border width in pixels (can be 1 or 0)
               GL_BGR,             // input image format
               imgDataType, // image data type (https://www.khronos.org/opengl/wiki/OpenGL_Type)
               tmp.ptr());  // pointer to data
  // glGenerateMipmap(GL_TEXTURE_2D);
}

void setWindowFullscreen() {
  // set window to fullscreen
  const ImGuiViewport *viewport = ImGui::GetMainViewport();
  ImGui::SetNextWindowPos(viewport->WorkPos);
  ImGui::SetNextWindowSize(viewport->WorkSize);

  ImGui::SetNextWindowViewport(viewport->ID);
  ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
  ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
}
