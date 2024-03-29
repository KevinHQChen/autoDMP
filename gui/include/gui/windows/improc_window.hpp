#pragma once

#include "gui/components/button.hpp"
#include "gui/components/slider.hpp"
// #include "gui/components/image.hpp"
#include "gui/guiframe.hpp"
#include "window.hpp"

#include "imcap/imcap.hpp"
#include "improc/improc.hpp"

namespace gui {

class ImProcWindow : public Window {
public:
  ImProcWindow(ImCap *imCap, ImProc *imProc);
  ~ImProcWindow();

  bool improcSetupVisible_{false}, improcVisible_{false}, startedImCap_{false},
      startedImProc_{false};

  void render() override;

private:
  void renderImCap();
  void renderImProc();
  void renderImProcConfigTable();

  // ImGuiWindowFlags imCapFlags = 0;
  // ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;

  ImCap *imCap_;
  ImProc *imProc_;

  // std::unique_ptr<IMMImage> rawImage_, procImage_;
  // std::array<std::unique_ptr<IMMImage>, numChans_> chImages_;
  std::unique_ptr<Toggle> imCapToggle_, imProcSetupToggle_, imProcToggle_;

  GUIFrame rawFrame;
  double scale{0.64}, scaleMin{0.1}, scaleMax{1.0};
  std::vector<GUIFrame> procGUIFrames;
  std::vector<double> y;

  std::vector<RotRect> chBBoxes;
  int chWidth, no_, bgSubHistory;
  double bgSubThres;

  cv::Point junction;
  ImVec2 mousePos, startPos, center, endPos, pointPos;
  ImDrawList *drawList;
  bool isDrawing{false}, drawChs{false}, drawJunc{false};

  ImVec2 clampToBounds(const ImVec2 &pos, float width, float height) {
    ImVec2 clamped;
    clamped.x = std::clamp(pos.x, 0.0f, width);
    clamped.y = std::clamp(pos.y, 0.0f, height);
    return clamped;
  }

  void draw(ImVec2 imgOrigin, double scale) {
    drawList = ImGui::GetWindowDrawList();
    // draw current junction position
    drawList->AddCircleFilled(
        ImVec2(imgOrigin.x + junction.x * scale, imgOrigin.y + junction.y * scale), 5.0f,
        IM_COL32(0, 255, 0, 255));
    // draw current chBBoxes and associated channel vectors
    for (const auto &chBBox : chBBoxes) {
      // draw current chBBoxes
      drawList->AddRect(ImVec2(imgOrigin.x + chBBox.x * scale, imgOrigin.y + chBBox.y * scale),
                        ImVec2(imgOrigin.x + (chBBox.x + chBBox.width) * scale,
                               imgOrigin.y + (chBBox.y + chBBox.height) * scale),
                        IM_COL32(0, 255, 0, 255));
      // draw associated channel vectors
      center = ImVec2(imgOrigin.x + (chBBox.x + chBBox.width / 2.0f) * scale,
                      imgOrigin.y + (chBBox.y + chBBox.height / 2.0f) * scale);
      float length = std::min(chBBox.width, chBBox.height) / 2.0f * scale;
      float radianAngle = (chBBox.angle * M_PI) / 180.0f;
      endPos = ImVec2(center.x - length * sin(radianAngle), center.y - length * cos(radianAngle));
      drawList->AddLine(center, endPos, IM_COL32(255, 255, 255, 255), 1.0f);
      drawList->AddCircleFilled(endPos, 3.0f, IM_COL32(255, 0, 0, 255));
    }

    // TODO make these 2 functions scale with the image (this is pretty annoying to do)
    drawNextChBBox(imgOrigin);
    drawNextJunc(imgOrigin);
  }

  void drawNextChBBox(ImVec2 imgOrigin) {
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

  void drawNextJunc(ImVec2 imgOrigin) {
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

  void drawFgLocs(int ch, ImVec2 imgOrigin, double y1, double y2) {
    drawList = ImGui::GetWindowDrawList();
    // draw primary and secondary fgLoc positions
    // info("drawing channel {} at {}, {}, offset by {}", ch, y1, y2, imgOrigin.y);
    drawList->AddCircleFilled(ImVec2(chWidth / 2.0 + imgOrigin.x, y1 + imgOrigin.y), 5.0f,
                              IM_COL32(0, 255, 0, 255));
    drawList->AddCircleFilled(ImVec2(chWidth / 2.0 + imgOrigin.x, y2 + imgOrigin.y), 5.0f,
                              IM_COL32(0, 255, 0, 255));
  }
};

} // namespace gui
