#include "gui/gui.hpp"

GUI::GUI()
    : conf(toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml")),
      guiConf(toml::find<guiConfig>(conf, "gui")), imCap(new ImCap()), imProc(new ImProc(imCap))
// , supervisor(new Supervisor())
{
  info("Config type: {}", type_name<decltype(guiConf)>());
  info("Parsed config: {}", toml::find(conf, "gui"));
}

GUI::~GUI() { delete imCap; }

void GUI::showRawImCap() {
  // set window to fullscreen
  const ImGuiViewport *viewport = ImGui::GetMainViewport();
  ImGui::SetNextWindowPos(viewport->WorkPos);
  ImGui::SetNextWindowSize(viewport->WorkSize);

  // show most recent frame or repeat previous frame if no new frame available
  rawFrame = imCap->getRawFrame();
  if (!rawFrame.empty()) {
    updateTexture(rawFrame, rawTextureID);
    rawWidth = rawFrame.cols;
    rawHeight = rawFrame.rows;
  }
  dear::Begin("Raw Image Capture", &guiConf.startImCap, imCapFlags) && [this]() {
    if (rawWidth > 0 && rawHeight > 0)
      ImGui::Image((void *)(intptr_t)rawTextureID, ImVec2(rawWidth, rawHeight));
    else
      ImGui::Text("Empty frame");
  };
}

void GUI::showTmplMatchSetup() {
  // TODO load template from file
  // if(toml::get<std::string>(conf["improc"]["tmplSrc"]) == "fromFile") {

  // }
  // // grab template from user-selected frame
  // else {
  //   cv::FileStorage chanPoseFile("chanPose.yml", cv::FileStorage::WRITE);
  //   preFrame = imCap->getPreFrame();
  //   if (!preFrame.empty()) {
  //     preWidth = preFrame.cols;
  //     preHeight = preFrame.rows;
  //   }
  //   info("Select bounding boxes for channel(s) (first channel must contain a droplet interface).");
  //   cv::selectROIs("channel selection", preFrame, chanPose.chanBBox, true, false);
  //   for (int i = 0; i < chanPose.chanBBox.size(); i++) {
  //     // save cropped channel as separate image
  //     currChan = preFrame(chanPose.chanBBox[i]).clone();
  //     // straighten any channels if necessary (for Y-junctions, the droplet will always start in an angled channel, thus we need to straighten the channel and define another bbox for it)

  //   }

  // }

  // imProc->setupTmplMatch();

  preFrame = imCap->getPreFrame();
  if (!preFrame.empty()) {
    updateTexture(preFrame, preTextureID);
    preWidth = preFrame.cols;
    preHeight = preFrame.rows;
  }

  ImGuiWindowFlags tmplMatchFlags = 0;
  tmplMatchFlags |= ImGuiWindowFlags_NoBackground;

  dear::Begin("Template Match Setup", &guiConf.setupTmplMatch, tmplMatchFlags) && [this]() {
    ImGui::Checkbox("Enable grid", &opt_enable_grid);
    ImGui::Checkbox("Enable context menu", &opt_enable_context_menu);
    ImGui::Checkbox("Draw rectangle", &opt_enable_rect);
    ImGui::Text(
        "Mouse Left: drag to add lines,\nMouse Right: drag to scroll, click for context menu.");

    // Using InvisibleButton() as a convenience 1) it will advance the layout cursor and 2) allows
    // us to use IsItemHovered()/IsItemActive()
    ImVec2 canvas_p0 = ImGui::GetCursorScreenPos();    // ImDrawList API uses screen coordinates!
    ImVec2 canvas_sz = ImGui::GetContentRegionAvail(); // Resize canvas to what's available
    if (canvas_sz.x < 50.0f)
      canvas_sz.x = 50.0f;
    if (canvas_sz.y < 50.0f)
      canvas_sz.y = 50.0f;
    ImVec2 canvas_p1 = ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y);

    // Draw border and background color
    ImGuiIO &io = ImGui::GetIO();
    ImDrawList *draw_list = ImGui::GetWindowDrawList();
    // draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));
    // Draw current frame
    // if (preWidth > 0 && preHeight > 0)
    //   ImGui::Image((void *)(intptr_t)preTextureID, ImVec2(preWidth, preHeight));
    if (preWidth > 0 && preHeight > 0)
      draw_list->AddImage((void *)(intptr_t)preTextureID, canvas_p0, canvas_p1);
    draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));

    // This will catch our interactions
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


    // Context menu (under default mouse threshold)
    ImVec2 drag_delta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Right);
    if (opt_enable_context_menu && drag_delta.x == 0.0f && drag_delta.y == 0.0f)
      ImGui::OpenPopupOnItemClick("context", ImGuiPopupFlags_MouseButtonRight);
    if (ImGui::BeginPopup("context"))
    {
      if (adding_line)
        points.resize(points.size() - 2);
      adding_line = false;
      if (ImGui::MenuItem("Remove one", NULL, false, points.Size > 0)) { points.resize(points.size() - 2); }
      if (ImGui::MenuItem("Remove all", NULL, false, points.Size > 0)) { points.clear(); }
      ImGui::EndPopup();
    }

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
    if (opt_enable_rect) {
      draw_list->AddLine(ImVec2(origin.x + rectStart.x, origin.y + rectStart.y),
                         ImVec2(origin.x + rectStart.x, origin.y + rectEnd.y),
                         IM_COL32(255, 255, 0, 255), 2.0f);
      draw_list->AddLine(ImVec2(origin.x + rectStart.x, origin.y + rectEnd.y),
                         ImVec2(origin.x + rectEnd.x, origin.y + rectEnd.y),
                         IM_COL32(255, 255, 0, 255), 2.0f);
      draw_list->AddLine(ImVec2(origin.x + rectEnd.x, origin.y + rectEnd.y),
                         ImVec2(origin.x + rectEnd.x, origin.y + rectStart.y),
                         IM_COL32(255, 255, 0, 255), 2.0f);
      draw_list->AddLine(ImVec2(origin.x + rectEnd.x, origin.y + rectStart.y),
                         ImVec2(origin.x + rectStart.x, origin.y + rectStart.y),
                         IM_COL32(255, 255, 0, 255), 2.0f);
    }
    draw_list->PopClipRect();
  };
}

void GUI::showImProc() {
  // show most recent frame or repeat previous frame if no new frame available
  procFrames = imProc->getProcFrames();
  int idx = 0;
  for (auto &frame : procFrames) {
    if (!frame.empty()) {
      updateTexture(frame, procTextureIDs[idx]);
      procWidths[idx] = frame.cols;
      procHeights[idx] = frame.rows;
    }
  }
  dear::Begin("Processed Image Capture", &guiConf.startImProc) && [this]() {
    int idx = 0;
    for (auto &textureID : procTextureIDs) {
      if (procWidths[idx] > 0 && procHeights[idx] > 0)
        ImGui::Image((void *)(intptr_t)textureID, ImVec2(procWidths[idx], procHeights[idx]));
      else
        ImGui::Text("Empty frame %d", idx);
    }
  };
}

ImGuiWrapperReturnType GUI::render() {
  dear::Begin("Menu") && [this]() {
    ImGui::Text("Instructions: TODO");
    ImGui::Text("Help (this might be a button?)");
    dear::MainMenuBar() && [this]() {
      dear::Menu("File") && [this]() { needToQuit = ImGui::MenuItem("Quit"); };
      dear::Menu("Setup") &&
          [this]() { ImGui::MenuItem("Image Capture", nullptr, &guiConf.startImCap);
        ImGui::MenuItem("Template Matching", nullptr, &guiConf.setupTmplMatch); };
      dear::Menu("Debug") &&
          [this]() { ImGui::MenuItem("Show Demo Window", nullptr, &guiConf.showDebug); };
    };
  };

  if (guiConf.startImCap) {
    if (!imCap->started())
      imCap->startCaptureThread();
    showRawImCap();
  } else {
    if (imCap->started())
      imCap->stopCaptureThread();
  }

  if (guiConf.setupTmplMatch) {
    if (imProc->started())
      imProc->stopImProcThread();
    showTmplMatchSetup();
  }

  // if (guiConf.startImProc) {
  //   if (!imProc->started())
  //     imProc->startImProcThread();
  //   showImProc();
  // } else {
  //   if (imProc->started())
  //     imProc->stopImProcThread();
  // }

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
