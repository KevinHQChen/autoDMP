#include "gui/windows/improc_window.hpp"

namespace gui {

ImProcWindow::ImProcWindow(ImCap *imCap, ImProc *imProc) : imCap_(imCap), imProc_(imProc) {
  imProcSetupToggle_ = std::make_unique<Toggle>("ImProc Setup", &improcSetupVisible_);
  // rawImage_ = std::make_unique<IMMImage>("Raw Image", 1);
  // procImage_ = std::make_unique<IMMImage>("Proc Image", 1);
  // for (int i = 0; i < NUM_TEMPLATES; i++)
  //   tmplImages_[i] = std::make_unique<IMMImage>("TmpImage: " + std::to_string(i), 1, true);
  // for (int i = 0; i < numChans_; i++)
  //   chImages_[i] = std::make_unique<IMMImage>("Channel " + std::to_string(i), 1);
}

ImProcWindow::~ImProcWindow() {
  imProcSetupToggle_.reset();
  // rawImage_.reset();
  // procImage_.reset();
  // for (int i = 0; i < NUM_TEMPLATES; i++)
  //   tmplImages_[i].reset();
  // for (int i = 0; i < numChans_; i++)
  //   chImages_[i].reset();
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
    numChs = chBBoxes.size();

    // Template BBoxes
    deleteFlags = std::vector<bool>(tmplBBoxes.size(), false);
    int prevRow = rowIndex;
    rowIndex = 0;
    for (auto &tmplBBox : tmplBBoxes) {
      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      std::string tmplLabel = "tmplBBox" + std::to_string(rowIndex);
      ImGui::Text("%s", tmplLabel.c_str());
      ImGui::TableSetColumnIndex(1);
      ImGui::PushID(prevRow + rowIndex * 5 + 0); // Unique ID for each editable field
      ImGui::InputInt("##x", &tmplBBox.x);
      ImGui::PopID();
      ImGui::TableSetColumnIndex(2);
      ImGui::PushID(prevRow + rowIndex * 5 + 1);
      ImGui::InputInt("##y", &tmplBBox.y);
      ImGui::PopID();
      ImGui::TableSetColumnIndex(3);
      ImGui::PushID(prevRow + rowIndex * 5 + 2);
      ImGui::InputInt("##width", &tmplBBox.width);
      ImGui::PopID();
      ImGui::TableSetColumnIndex(4);
      ImGui::PushID(prevRow + rowIndex * 5 + 3);
      ImGui::InputInt("##height", &tmplBBox.height);
      ImGui::PopID();
      ImGui::TableSetColumnIndex(6);
      std::string deleteBtnLbl = "Delete Tmpl " + std::to_string(rowIndex);
      if (ImGui::Button(deleteBtnLbl.c_str()))
        deleteFlags[rowIndex] = true;
      rowIndex++;
    }
    // Delete flagged items
    for (int i = deleteFlags.size() - 1; i >= 0; --i)
      if (deleteFlags[i])
        tmplBBoxes.erase(tmplBBoxes.begin() + i);
    numTmpls = 2 * tmplBBoxes.size();

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
  ImGui::InputInt("Num Templates: ", &numTmpls);
  ImGui::InputInt("Num Channels: ", &numChs);
  ImGui::SliderInt("Channel Width: ", &chWidth, 0, 100);
  ImGui::SliderFloat("Template Threshold: ", &tmplThres, 0.0f, 1.0f, "ratio = %.3f");

  ImGui::Text("Image dimensions: %dx%d", rawFrame.width, rawFrame.height);
}

void ImProcWindow::render() {
  if (ImGui::IsKeyPressed(ImGuiKey_C))
    !imCap_->started() ? imCap_->startCaptureThread() : imCap_->stopCaptureThread();
  if (ImGui::IsKeyPressed(ImGuiKey_Enter))
    visible_ = !visible_;
  if (visible_) {
    if (ImGui::Begin("Raw Image Capture", &visible_)) {
      // draw image
      rawFrame = imCap_->getFrame();
      ImGui::Image((ImTextureID)rawFrame.texture, ImVec2(rawFrame.width, rawFrame.height));
      imgOrigin = ImGui::GetItemRectMin();

      // improc setup
      if (ImGui::IsKeyPressed(ImGuiKey_S))
        improcSetupVisible_ = !improcSetupVisible_;
      if (improcSetupVisible_) {
        if (ImGui::Button("Select Junction"))
          drawJunc = true;
        ImGui::SameLine();
        if (ImGui::Button("Select Channel"))
          drawChs = true;
        ImGui::SameLine();
        if (ImGui::Button("Select Template"))
          drawTmpl = true;
        draw();
        renderImProcConfigTable();

        if (ImGui::Button("Update ImProc Config")) {
          imProc_->impConf.setChROIs(chBBoxes);
          imProc_->impConf.chWidth_ = chWidth;
          imProc_->impConf.numChs_ = numChs;
          imProc_->impConf.tmplThres_ = tmplThres;
          imProc_->impConf.numTmpls_ = numTmpls;
          imProc_->impConf.clearTmplImgs();
          for (const auto &tmplBBox : tmplBBoxes) {
            cv::Mat tmplImg, tmplImgFlip;
            tmplImg = rawFrame.mat(tmplBBox).clone();
            imProc_->impConf.setTmplImg(tmplImg);
            cv::flip(tmplImg, tmplImgFlip, -1); // 180deg CCW (flip around x & y-axis)
            imProc_->impConf.setTmplImg(tmplImgFlip);
          }
        }
        ImGui::SameLine();
        if (ImGui::Button("Save ImProc Config"))
          imProc_->saveConfig();
        ImGui::SameLine();
        if (ImGui::Button("Load ImProc Config")) {
          imProc_->loadConfig();
          chBBoxes = imProc_->impConf.getChROIs();
          chWidth = imProc_->impConf.chWidth_;
          numChs = imProc_->impConf.numChs_;
          tmplThres = imProc_->impConf.tmplThres_;
          numTmpls = imProc_->impConf.numTmpls_;
          tmplBBoxes.clear();
        }

      }
      ImGui::End();
    }
  }

  if (ImGui::IsKeyPressed(ImGuiKey_I))
    improcVisible_ = !improcVisible_;
  if (improcVisible_) {
    imProc_->startProcThread();
    if (ImGui::Begin("Channels", &improcVisible_)) {
      for (int i = 0; i < imProc_->impConf.numChs_; ++i) {
        if (i != 0)
          ImGui::SameLine();
        procGUIFrames[i] = imProc_->getProcFrame(i);
        ImGui::Image((ImTextureID)procGUIFrames[i].texture,
                     ImVec2(procGUIFrames[i].width, procGUIFrames[i].height));
      }
      ImGui::End();
    }
    if (ImGui::Begin("Proc Image", &improcVisible_)) {
      preFrame = imProc_->getProcFrame();
      ImGui::Image((ImTextureID)preFrame.texture, ImVec2(preFrame.width, preFrame.height));
      ImGui::End();
    }
  } else
    imProc_->stopProcThread();
}

} // namespace gui
