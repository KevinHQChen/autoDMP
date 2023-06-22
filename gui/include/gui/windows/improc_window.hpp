#pragma once

#include "gui/components/button.hpp"
#include "gui/components/slider.hpp"
// #include "gui/components/image.hpp"
#include "gui/guiframe.hpp"
#include "window.hpp"

#include "imcap/imcap.hpp"
#include "improc/improc.hpp"

namespace gui {

inline ImVec2 clampToBounds(const ImVec2 &pos, float width, float height) {
  ImVec2 clamped;
  clamped.x = std::clamp(pos.x, 0.0f, width);
  clamped.y = std::clamp(pos.y, 0.0f, height);
  return clamped;
}

class ImProcWindow : public Window {
  // ImGuiWindowFlags imCapFlags = 0;
  // ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;

  ImCap *imCap_;
  ImProc *imProc_;

  std::unique_ptr<Toggle> imProcSetupToggle_;

  // std::unique_ptr<IMMImage> rawImage_, procImage_;
  // std::array<std::unique_ptr<IMMImage>, numChans_> chImages_;

  GUIFrame rawFrame, preFrame;
  GUIFrame procGUIFrames[3];

  std::vector<RotRect> chBBoxes;
  int chWidth, numChs, bgSubHistory;
  double bgSubThres;

  cv::Point junction;
  ImVec2 mousePos, startPos, center, endPos, pointPos, imgOrigin;
  ImDrawList *drawList;
  bool isDrawing{false}, drawChs{false}, drawJunc{false};

  void draw() {
    drawList = ImGui::GetWindowDrawList();
    // draw current junction position
    drawList->AddCircleFilled(ImVec2(junction.x + imgOrigin.x, junction.y + imgOrigin.y), 5.0f,
                              IM_COL32(0, 255, 0, 255));
    // draw current chBBoxes and associated channel vectors
    for (const auto &chBBox : chBBoxes) {
      // draw current chBBoxes
      drawList->AddRect(
          ImVec2(chBBox.x + imgOrigin.x, chBBox.y + imgOrigin.y),
          ImVec2(chBBox.x + chBBox.width + imgOrigin.x, chBBox.y + chBBox.height + imgOrigin.y),
          IM_COL32(0, 255, 0, 255));
      // draw associated channel vectors
      center = ImVec2(chBBox.x + chBBox.width / 2.0f + imgOrigin.x,
                      chBBox.y + chBBox.height / 2.0f + imgOrigin.y);
      float length = std::min(chBBox.width, chBBox.height) / 2.0f;
      float radianAngle = (chBBox.angle * M_PI) / 180.0f;
      endPos = ImVec2(center.x - length * sin(radianAngle), center.y - length * cos(radianAngle));
      drawList->AddLine(center, endPos, IM_COL32(255, 255, 255, 255), 1.0f);
      drawList->AddCircleFilled(endPos, 3.0f, IM_COL32(255, 0, 0, 255));
    }

    drawNextChBBox();
    drawNextJunc();
  }

  void drawNextChBBox() {
    if (drawChs) {
      mousePos = ImGui::GetIO().MousePos;
      if (ImGui::IsMouseDown(0) && !isDrawing) {
        startPos = ImVec2(mousePos.x - imgOrigin.x, mousePos.y - imgOrigin.y);
        startPos = clampToBounds(startPos, rawFrame.width, rawFrame.height);
        isDrawing = true;
      }
      if (isDrawing) {
        endPos = ImVec2(mousePos.x - imgOrigin.x, mousePos.y - imgOrigin.y);
        endPos = clampToBounds(endPos, rawFrame.width, rawFrame.height);
        drawList->AddRect(ImVec2(startPos.x + imgOrigin.x, startPos.y + imgOrigin.y),
                          ImVec2(endPos.x + imgOrigin.x, endPos.y + imgOrigin.y),
                          IM_COL32(255, 0, 0, 255));
      }
      if (ImGui::IsMouseReleased(0) && isDrawing) {
        endPos = ImVec2(mousePos.x - imgOrigin.x, mousePos.y - imgOrigin.y);
        endPos = clampToBounds(endPos, rawFrame.width, rawFrame.height);
        chBBoxes.push_back(cv::Rect((startPos.x < endPos.x) ? startPos.x : endPos.x,
                                    (startPos.y < endPos.y) ? startPos.y : endPos.y,
                                    std::abs(endPos.x - startPos.x),
                                    std::abs(endPos.y - startPos.y)));
        isDrawing = drawChs = false;
      }
    }
  }

  void drawNextJunc() {
    if (drawJunc) {
      mousePos = ImGui::GetIO().MousePos;
      if (ImGui::IsMouseClicked(0)) {
        pointPos = ImVec2(mousePos.x - imgOrigin.x, mousePos.y - imgOrigin.y);
        pointPos = clampToBounds(pointPos, rawFrame.width, rawFrame.height);
        junction = cv::Point(pointPos.x, pointPos.y);
        drawJunc = false;
      }
    }
  }

  void renderImCap();
  void renderImProc();
  void renderImProcConfigTable();

public:
  bool improcSetupVisible_{false}, improcVisible_{false};

  ImProcWindow(ImCap *imCap, ImProc *imProc);
  ~ImProcWindow();
  void render() override;
};

} // namespace gui
