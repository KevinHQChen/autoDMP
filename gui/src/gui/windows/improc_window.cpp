#include "gui/windows/improc_window.hpp"

namespace gui {

ImProcWindow::ImProcWindow(ImCap *imCap, ImProc *imProc) : imCap_(imCap), imProc_(imProc) {
  // rawImage_ = std::make_unique<IMMImage>("Raw Image", 1);
  // procImage_ = std::make_unique<IMMImage>("Proc Image", 1);
  // for (int i = 0; i < numChans_; i++)
  //   chImages_[i] = std::make_unique<IMMImage>("Channel " + std::to_string(i), 1);
  imCapToggle_ = std::make_unique<Toggle>("ImCap (c)", &startedImCap_);
  imProcSetupToggle_ = std::make_unique<Toggle>("ImProc Setup (s)", &improcSetupVisible_);
  imProcToggle_ = std::make_unique<Toggle>("ImProc (i)", &startedImProc_);
}

ImProcWindow::~ImProcWindow() {
  // rawImage_.reset();
  // procImage_.reset();
  // for (int i = 0; i < numChans_; i++)
  //   chImages_[i].reset();
  imCapToggle_.reset();
  imProcSetupToggle_.reset();
  imProcToggle_.reset();
}

void ImProcWindow::render() {
  renderImCap();
  renderImProc();
}

void ImProcWindow::renderImCap() {
  if (ImGui::IsKeyPressed(ImGuiKey_C) || imCapToggle_->changed())
    !imCap_->started() ? imCap_->startThread() : imCap_->stopThread();
  if (ImGui::IsKeyPressed(ImGuiKey_D))
    visible_ = !visible_;
  if (ImGui::IsKeyPressed(ImGuiKey_S))
    improcSetupVisible_ = !improcSetupVisible_;
  if (visible_ && ImGui::Begin("Raw Image Capture", &visible_)) {
    // draw image
    if (imCap_->started()) {
      startedImCap_ = true;
      rawFrame = imCap_->getFrame();
    } else
      startedImCap_ = false;
    ImGui::Image((ImTextureID)rawFrame.texture, ImVec2(rawFrame.width, rawFrame.height));
    ImVec2 imgOrigin = ImGui::GetItemRectMin();

    // improc setup
    imCapToggle_->render();
    imProcSetupToggle_->render();
    if (improcSetupVisible_) {
      if (ImGui::Button("Select Junction"))
        drawJunc = true;
      ImGui::SameLine();
      if (ImGui::Button("Select Channel"))
        drawChs = true;
      draw(imgOrigin);
      renderImProcConfigTable();

      if (ImGui::Button("Update ImProc Config")) {
        imProc_->impConf.setChROIs(chBBoxes);
        imProc_->impConf.setChWidth(chWidth);
        imProc_->impConf.setNumChs(no_);
        imProc_->impConf.setBgSubHistory(bgSubHistory);
        imProc_->impConf.setBgSubThres(bgSubThres);
        procGUIFrames.clear();
        for (int i = 0; i < no_; ++i) {
          GUIFrame emptyFrame;
          emptyFrame = cv::Mat(0, 0, CV_8UC1);
          procGUIFrames.push_back(emptyFrame);
        }
      }
      ImGui::SameLine();
      if (ImGui::Button("Save ImProc Config"))
        imProc_->saveConfig();
      ImGui::SameLine();
      if (ImGui::Button("Load ImProc Config")) {
        imProc_->loadConfig();
        chBBoxes = imProc_->impConf.getChROIs();
        chWidth = imProc_->impConf.getChWidth();
        no_ = imProc_->impConf.getNumChs();
        bgSubHistory = imProc_->impConf.getBgSubHistory();
        bgSubThres = imProc_->impConf.getBgSubThres();
      }
    }
    ImGui::End();
  }
}

void ImProcWindow::renderImProc() {
  if (ImGui::IsKeyPressed(ImGuiKey_I) || imProcToggle_->changed())
    !imProc_->started() ? imProc_->startThread() : imProc_->stopThread();
  if (ImGui::IsKeyPressed(ImGuiKey_J))
    improcVisible_ = !improcVisible_;
  if (improcVisible_ && ImGui::Begin("Channels", &improcVisible_)) {
    y = imProc_->getY();
    if (!y.empty()) {
      for (int i = 0; i < no_; ++i) {
        if (i != 0)
          ImGui::SameLine();
        if (imProc_->started()) {
          startedImProc_ = true;
          procGUIFrames[i] = imProc_->getProcFrame(i);
        } else
          startedImProc_ = false;
        ImGui::Image((ImTextureID)procGUIFrames[i].texture,
                     ImVec2(procGUIFrames[i].width, procGUIFrames[i].height));
        drawFgLocs(i, ImGui::GetItemRectMin(), -y[i] + imProc_->impConf.getChWidth() / 2.0,
                   -y[i + no_] + imProc_->impConf.getChWidth() / 2.0);
        // print y1 and y2
        ImGui::SameLine();
        ImGui::Text("y1: %.1f", y[i]);
        ImGui::SameLine();
        ImGui::Text("y2: %.1f", y[i + no_]);
      }
    }
    imProcToggle_->render();
    ImGui::End();
  }
}

void ImProcWindow::renderImProcConfigTable() {
  // Add table to display and edit properties
  if (ImGui::BeginTable("prop_table", 7)) {
    ImGui::TableSetupColumn("Property");
    ImGui::TableSetupColumn("X");
    ImGui::TableSetupColumn("Y");
    ImGui::TableSetupColumn("Width");
    ImGui::TableSetupColumn("Height");
    ImGui::TableSetupColumn("Rot Angle");
    ImGui::TableSetupColumn("Delete");
    ImGui::TableHeadersRow();

    // Channel BBoxes
    std::vector<bool> deleteFlags(chBBoxes.size(), false);
    int rowIndex = 0;
    for (auto &chBBox : chBBoxes) {
      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      std::string chLabel = "chBBox" + std::to_string(rowIndex);
      ImGui::Text("%s", chLabel.c_str());
      ImGui::TableSetColumnIndex(1);
      ImGui::PushID(rowIndex * 5 + 0); // Unique ID for each editable field
      ImGui::InputInt("##x", &chBBox.x);
      ImGui::PopID();
      ImGui::TableSetColumnIndex(2);
      ImGui::PushID(rowIndex * 5 + 1);
      ImGui::InputInt("##y", &chBBox.y);
      ImGui::PopID();
      ImGui::TableSetColumnIndex(3);
      ImGui::PushID(rowIndex * 5 + 2);
      ImGui::InputInt("##width", &chBBox.width);
      ImGui::PopID();
      ImGui::TableSetColumnIndex(4);
      ImGui::PushID(rowIndex * 5 + 3);
      ImGui::InputInt("##height", &chBBox.height);
      ImGui::PopID();
      ImGui::TableSetColumnIndex(5);
      ImGui::PushID(rowIndex * 5 + 4);
      ImGui::InputInt("##rotAngle", &chBBox.angle);
      ImGui::PopID();
      ImGui::TableSetColumnIndex(6);
      std::string deleteBtnLbl = "Delete Ch " + std::to_string(rowIndex);
      if (ImGui::Button(deleteBtnLbl.c_str()))
        deleteFlags[rowIndex] = true;
      rowIndex++;
    }
    // Delete flagged items
    for (int i = deleteFlags.size() - 1; i >= 0; --i)
      if (deleteFlags[i])
        chBBoxes.erase(chBBoxes.begin() + i);
    no_ = chBBoxes.size();

    // Junction
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::Text("Junction");
    ImGui::TableSetColumnIndex(1);
    ImGui::InputInt("##Junction_x", &junction.x);
    ImGui::TableSetColumnIndex(2);
    ImGui::InputInt("##Junction_y", &junction.y);
    ImGui::EndTable();
  }
  ImGui::InputInt("Channel Width: ", &chWidth);
  ImGui::InputInt("Num Channels: ", &no_);
  ImGui::InputInt("BgSub History", &bgSubHistory);
  ImGui::InputDouble("BgSub Threshold", &bgSubThres);

  ImGui::Text("Image dimensions: %dx%d", rawFrame.width, rawFrame.height);
}

} // namespace gui
