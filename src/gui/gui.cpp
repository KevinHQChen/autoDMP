#include "gui/gui.hpp"

GUI::GUI()
    : conf(TOML11_PARSE_IN_ORDER("config/setup.toml")), guiConf(toml::find<guiConfig>(conf, "gui")),
      imCap(new ImCap()), imProc(new ImProc(imCap)), supervisor(new Supervisor(imProc)) {
  info("Config type: {}", type_name<decltype(guiConf)>());
  info("Parsed config: {}", toml::find(conf, "gui"));
  // TODO may want to make GUI, ImCap, ImProc, Supervisor, etc. singletons
}

GUI::~GUI() {
  delete supervisor;
  delete imProc;
  delete imCap;
}

void GUI::showRawImCap() {
  if (guiConf.startImCap) {
    imCap->startCaptureThread();

    if (ImGui::Begin("Raw Image Capture", &guiConf.startImCap, imCapFlags)) {
      rawFrame = imCap->getRawFrame();
      (rawFrame.empty) ? ImGui::Text("No image available")
                       : ImGui::Image((ImTextureID)rawFrame.texture,
                                      ImVec2(rawFrame.width, rawFrame.height));
      ImGui::End();
    }
  } else
    imCap->stopCaptureThread();
}

void GUI::showImProc() {
  if (guiConf.startImProc) {
    imProc->startProcThread();
    for (int idx = 0; idx < toml::get<int>(conf["improc"]["numChans"]); ++idx) {
      std::string windowName = "Processed Frame " + std::to_string(idx);
      if (ImGui::Begin(windowName.c_str(), &guiConf.startImProc)) {
        procGUIFrames[idx] = imProc->getProcFrame(idx);
        (procGUIFrames[idx].empty)
            ? ImGui::Text("Empty frame %d", idx)
            : ImGui::Image((ImTextureID)procGUIFrames[idx].texture,
                           ImVec2(procGUIFrames[idx].width, procGUIFrames[idx].height));
        ImGui::End();
      }
    }

    if (ImGui::Begin("Proc Image Capture", &guiConf.startImProc)) {
      preFrame = imProc->getProcFrame();
      (preFrame.empty) ? ImGui::Text("No image available")
                       : ImGui::Image((ImTextureID)preFrame.texture,
                                      ImVec2(preFrame.width, preFrame.height));
      ImGui::End();
    }
  } else
    imProc->stopProcThread();
}

void GUI::showImProcSetup() {
  imProc->setSetupStatus(guiConf.startImProcSetup);
  if (guiConf.startImProcSetup) {
    if (ImGui::Begin("ImProc Setup", &guiConf.startImProcSetup)) {
      // update bbox
      int bbox[4] = {imProc->impConf.getBBox().x, imProc->impConf.getBBox().y,
                     imProc->impConf.getBBox().width,
                     imProc->impConf.getBBox().height}; // x, y, width, height
      ImGui::SliderInt("BBox.x", &bbox[0], 0, 1000);
      ImGui::SliderInt("BBox.y", &bbox[1], 0, 1000);
      ImGui::SliderInt("BBox.width", &bbox[2], 0, 1000);
      ImGui::SliderInt("BBox.height", &bbox[3], 0, 1000);
      imProc->impConf.setBBox(cv::Rect(bbox[0], bbox[1], bbox[2], bbox[3]));

      // update junction
      int junc[2] = {imProc->impConf.getJunction().x, imProc->impConf.getJunction().y};
      junc[0] = bbox[2] / 2;
      ImGui::SliderInt("junction.y", &junc[1], 0, 1000);
      imProc->impConf.setJunction(cv::Point(junc[0], junc[1]));

      // update chanWidth
      int chanWidth = imProc->impConf.getChanWidth();
      ImGui::SliderInt("Channel Width", &chanWidth, 0, 100);
      imProc->impConf.setChanWidth(chanWidth);

      // update rotAngles
      std::vector<int> rotAngles = imProc->impConf.getRotAngle();
      for (int idx = 0; idx < toml::find<int>(conf["improc"], "numChans"); ++idx) {
        std::string chanWinName = "Channel " + std::to_string(idx);
        if (ImGui::CollapsingHeader(chanWinName.c_str())) {
          std::string rotWinName = "Rot Angle " + std::to_string(idx);
          ImGui::SliderInt(rotWinName.c_str(), &rotAngles[idx], -180, 180);
        }
      }
      imProc->impConf.setRotAngle(rotAngles);

      // update chanBBox, rotChanBBox (using bbox, junction, chanWidth, rotAngle)
      std::vector<cv::Rect> chanBBoxes = imProc->impConf.getChanBBox();
      chanBBoxes[0] = cv::Rect(junc[0], junc[1], bbox[2] / 2, bbox[2] / 2);
      chanBBoxes[1] = cv::Rect(0, junc[1], bbox[2] / 2, bbox[2] / 2);
      chanBBoxes[2] = cv::Rect(junc[0] - chanWidth / 2, 0, chanWidth, junc[1]);
      imProc->impConf.setChanBBox(chanBBoxes);
      std::vector<cv::Rect> rotChanBBoxes = imProc->impConf.getRotChanBBox();
      rotChanBBoxes[0] = cv::Rect(bbox[2] / 2.0 * 1.414 / 2.0 - chanWidth / 2.0, 0, chanWidth,
                                  bbox[2] / 2.0 * 1.414);
      rotChanBBoxes[1] = cv::Rect(bbox[2] / 2.0 * 1.414 / 2.0 - chanWidth / 2.0, 0, chanWidth,
                                  bbox[2] / 2.0 * 1.414);
      rotChanBBoxes[2] = cv::Rect(0, 0, 0, 0);
      imProc->impConf.setRotChanBBox(rotChanBBoxes);

      // update tmplBBox
      int tmplBBox[4] = {imProc->impConf.getTmplBBox().x, imProc->impConf.getTmplBBox().y,
                         imProc->impConf.getTmplBBox().width,
                         imProc->impConf.getTmplBBox().height}; // x, y, width, height
      ImGui::SliderInt("TmplBBox.x", &tmplBBox[0], 0, 100);
      ImGui::SliderInt("TmplBBox.y", &tmplBBox[1], 0, 1000);
      ImGui::SliderInt("TmplBBox.width", &tmplBBox[2], 0, 100);
      ImGui::SliderInt("TmplBBox.height", &tmplBBox[3], 0, 100);
      imProc->impConf.setTmplBBox(cv::Rect(tmplBBox[0], tmplBBox[1], tmplBBox[2], tmplBBox[3]));

      // update tmplThres
      float tmplThres = imProc->tmplThres;
      ImGui::SliderFloat("Tmpl Thres", &tmplThres, 0.0f, 1.0f, "ratio = %.3f");
      imProc->tmplThres = tmplThres;

      // show tmplImg
      std::array<cv::Mat, NUM_TEMPLATES> tmplImg = imProc->impConf.getTmplImg();
      for (int idx = 0; idx < NUM_TEMPLATES; ++idx) {
        tmplGUIFrames[idx] = tmplImg[idx];
        (tmplGUIFrames[idx].empty)
            ? ImGui::Text("Empty tmpl %d", idx)
            : ImGui::Image((void *)(intptr_t)tmplGUIFrames[idx].texture,
                           ImVec2(tmplGUIFrames[idx].width, tmplGUIFrames[idx].height));
      }

      // load from/save to file
      if (ImGui::Button("Load from file"))
        imProc->loadConfig();
      if (ImGui::Button("Save to file"))
        imProc->saveConfig();
      if (ImGui::Button("Quit ImProc Setup"))
        guiConf.startImProcSetup = false;

      ImGui::End();
    }
  }
}

std::optional<int> GUI::render() {
  setWindowFullscreen();

  ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
  ImGui::Begin("DockSpace", nullptr, dockSpaceFlags); // not guaranteed to return true
  ImGui::PopStyleVar();
  ImGui::PopStyleVar(2);
  // Submit the DockSpace
  ImGuiID dockspace_id = ImGui::GetID("MyDockSpace");
  ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockNodeFlags);
  // Add menu to DockSpace
  if (ImGui::BeginMenuBar()) {
    if (ImGui::BeginMenu("File")) {
      needToQuit = ImGui::MenuItem("Quit");
      ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Setup")) {
      ImGui::MenuItem("Start Image Capture", nullptr, &guiConf.startImCap);
      ImGui::MenuItem("Setup Image Processing", nullptr, &guiConf.startImProcSetup);
      ImGui::MenuItem("Start Image Processing", nullptr, &guiConf.startImProc);
      ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Debug")) {
      ImGui::MenuItem("Show Demo Window", nullptr, &guiConf.showDebug);
      ImGui::EndMenu();
    }
    ImGui::EndMenuBar();
  }
  ImGui::End();

  showRawImCap();
  showImProcSetup();
  showImProc();
  if (guiConf.showDebug)
    ImGui::ShowDemoWindow(&guiConf.showDebug);
  if (needToQuit)
    return 0;

  return {};
}

void GUI::startGUIThread() {
  guiThread = std::thread(&GUI::imguiMain, this);
  guiThread.join();
}

static void glfw_error_callback(int error, const char *description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

/*
imgui_main initializes an ImGui openGL/glfw backend and then runs the passed
std::function<std::optional<int>()> is called repeatedly until the std::optional it returns has a
value, which is then returned as the exit code.
*/
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

void GUI::contextMenu(bool enable) {
  // Context menu (under default mouse threshold)
  ImVec2 drag_delta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Right);
  if (enable && drag_delta.x == 0.0f && drag_delta.y == 0.0f)
    ImGui::OpenPopupOnItemClick("context", ImGuiPopupFlags_MouseButtonRight);
  if (ImGui::BeginPopup("context")) {
    if (adding_line)
      points.resize(points.size() - 2);
    adding_line = false;
    if (ImGui::MenuItem("Remove one", NULL, false, points.Size > 0)) {
      points.resize(points.size() - 2);
    }
    if (ImGui::MenuItem("Remove all", NULL, false, points.Size > 0)) {
      points.clear();
    }
    ImGui::EndPopup();
  }
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

/*
void GUI::showImProcSetup() {
  // TODO load template from file
  if(toml::get<std::string>(conf["improc"]["tmplSrc"]) == "fromFile") {
    imProc->loadConfig(toml::get<std::string>(conf["improc"]["path"]));

  }
  // // grab template from user-selected frame
  // else {
  //   cv::FileStorage chanPoseFile("chanPose.yml", cv::FileStorage::WRITE);
  //   preFrame = imCap->getPreFrame();
  //   if (!preFrame.empty()) {
  //     preWidth = preFrame.cols;
  //     preHeight = preFrame.rows;
  //   }
  //   info("Select bounding boxes for channel(s) (first channel must contain a droplet
  //   interface)."); cv::selectROIs("channel selection", preFrame, chanPose.chanBBox, true, false);
  //   for (int i = 0; i < chanPose.chanBBox.size(); i++) {
  //     // save cropped channel as separate image
  //     currChan = preFrame(chanPose.chanBBox[i]).clone();
  //     // straighten any channels if necessary (for Y-junctions, the droplet will always start in
  //     an angled channel, thus we need to straighten the channel and define another bbox for it)

  //   }

  // }

  // imProc->setupTmplMatch();

  if (guiConf.setupImProc) {
    imProc->startSetupThread();

    if (ImGui::Begin("Image Processing Setup", &guiConf.setupImProc, imProcSetupFlags)) {
      ImGui::Checkbox("Enable grid", &opt_enable_grid);
      ImGui::Checkbox("Enable context menu", &opt_enable_context_menu);
      ImGui::Checkbox("Draw rectangle", &opt_enable_rect);
      ImGui::Text(
          "Mouse Left: drag to add lines,\nMouse Right: drag to scroll, click for context menu.");

      canvas_p0 = ImGui::GetCursorScreenPos();    // ImDrawList API uses screen coordinates!
      canvas_sz = ImGui::GetContentRegionAvail(); // Resize canvas to what's available
      if (canvas_sz.x < 50.0f)
        canvas_sz.x = 50.0f;
      if (canvas_sz.y < 50.0f)
        canvas_sz.y = 50.0f;
      canvas_p1 = ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y);

      // Draw border and background color
      ImGuiIO &io = ImGui::GetIO();
      ImDrawList *draw_list = ImGui::GetWindowDrawList();
      // draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));
      // Draw current frame
      preFrame = imCap->getPreFrame();
      (preFrame.empty)
          ? draw_list->AddText(canvas_p0, IM_COL32(255, 255, 255, 255), "No image available")
          : draw_list->AddImage((void *)(intptr_t)preFrame.texture, canvas_p0, canvas_p1);

      draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));

      // This will catch our interactions
      // Using InvisibleButton() as a convenience 1) it will advance the layout cursor and 2) allows
      // us to use IsItemHovered()/IsItemActive()
      ImGui::InvisibleButton("canvas", canvas_sz,
                             ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
      const bool is_hovered = ImGui::IsItemHovered(); // Hovered
      const bool is_active = ImGui::IsItemActive();   // Held
      const ImVec2 origin(canvas_p0.x + scrolling.x,
                          canvas_p0.y + scrolling.y); // Lock scrolled origin
      const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

      // Add first and second point
      if (is_hovered && !adding_line && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
        points.push_back(mouse_pos_in_canvas);
        points.push_back(mouse_pos_in_canvas);
        adding_line = true;
      }
      if (adding_line) {
        points.back() = mouse_pos_in_canvas;
        if (!ImGui::IsMouseDown(ImGuiMouseButton_Left))
          adding_line = false;
      }

      // Add rectangle
      if (opt_enable_rect) {
        if (is_hovered && !addingRect && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
          rectStart = mouse_pos_in_canvas;
          addingRect = true;
        }
        if (addingRect) {
          rectEnd = mouse_pos_in_canvas;
          if (!ImGui::IsMouseDown(ImGuiMouseButton_Left))
            addingRect = false;
        }
      }

      contextMenu(opt_enable_context_menu);

      // Draw grid, all lines, rectangle in the canvas
      draw_list->PushClipRect(canvas_p0, canvas_p1, true);
      // Draw grid
      if (opt_enable_grid) {
        const float GRID_STEP = 64.0f;
        for (float x = fmodf(scrolling.x, GRID_STEP); x < canvas_sz.x; x += GRID_STEP)
          draw_list->AddLine(ImVec2(canvas_p0.x + x, canvas_p0.y),
                             ImVec2(canvas_p0.x + x, canvas_p1.y), IM_COL32(200, 200, 200, 40));
        for (float y = fmodf(scrolling.y, GRID_STEP); y < canvas_sz.y; y += GRID_STEP)
          draw_list->AddLine(ImVec2(canvas_p0.x, canvas_p0.y + y),
                             ImVec2(canvas_p1.x, canvas_p0.y + y), IM_COL32(200, 200, 200, 40));
      }
      // Draw all lines
      for (int n = 0; n < points.Size; n += 2)
        draw_list->AddLine(ImVec2(origin.x + points[n].x, origin.y + points[n].y),
                           ImVec2(origin.x + points[n + 1].x, origin.y + points[n + 1].y),
                           IM_COL32(255, 255, 0, 255), 2.0f);
      // Draw rectangle
      if (opt_enable_rect)
        draw_list->AddRect(ImVec2(origin.x + rectStart.x, origin.y + rectStart.y),
                           ImVec2(origin.x + rectEnd.x, origin.y + rectEnd.y),
                           IM_COL32(255, 255, 0, 255), 0.0f, 0, 2.0f);
      draw_list->PopClipRect();
      ImGui::End();
    }
  } else
    imProc->stopSetupThread();
}
*/
