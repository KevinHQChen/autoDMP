#pragma once

#include "imgui.h"
#include "implot/implot.h"
#include "util/util.hpp"
#include <vector>

// utility structure for realtime plot
struct ScrollingBuffer {
  int MaxSize;
  int Offset;
  ImVector<ImVec2> Data;
  ScrollingBuffer(int max_size = 4000) {
    MaxSize = max_size;
    Offset = 0;
    Data.reserve(MaxSize);
  }
  void AddPoint(float x, float y) {
    if (Data.size() < MaxSize)
      Data.push_back(ImVec2(x, y));
    else {
      Data[Offset] = ImVec2(x, y);
      Offset = (Offset + 1) % MaxSize;
    }
  }
  void Erase() {
    if (Data.size() > 0) {
      Data.shrink(0);
      Offset = 0;
    }
  }
};

inline void plotVectorNN(std::string plotName, const char *xAx, const char *yAx, double yMin,
                         double yMax, std::vector<std::string> dataNames,
                         const std::vector<std::vector<ScrollingBuffer *>> &vecs, float guiTime,
                         float history) {
  if (ImPlot::BeginPlot(plotName.c_str(), ImVec2(-1, 300))) {
    ImPlot::SetupAxes(xAx, yAx); //, implotFlags, implotFlags);
    ImPlot::SetupAxisLimits(ImAxis_X1, guiTime - history, guiTime, ImGuiCond_Always);
    ImPlot::SetupAxisLimits(ImAxis_Y1, yMin, yMax);
    std::string label;
    for (int i = 0; i < dataNames.size(); ++i) {
      for (int j = 0; j < vecs[i].size(); ++j) {
        label = dataNames[i] + std::to_string(j);
        ImPlot::PlotLine(label.c_str(), &vecs[i][j]->Data[0].x, &vecs[i][j]->Data[0].y,
                         vecs[i][j]->Data.size(), 0, vecs[i][j]->Offset, 2 * sizeof(float));
      }
    }
    ImPlot::EndPlot();
  }
}

inline void plotVectorN(std::string plotName, const char *xAx, const char *yAx, double yMin,
                        double yMax, std::string dataName,
                        const std::vector<ScrollingBuffer *> &vecs, float guiTime, float history) {
  if (ImPlot::BeginPlot(plotName.c_str(), ImVec2(-1, 300))) {
    int idx = 0;
    ImPlot::SetupAxes(xAx, yAx); //, implotFlags, implotFlags);
    ImPlot::SetupAxisLimits(ImAxis_X1, guiTime - history, guiTime, ImGuiCond_Always);
    ImPlot::SetupAxisLimits(ImAxis_Y1, yMin, yMax);
    for (auto &vec : vecs)
      ImPlot::PlotLine((dataName + std::to_string(idx++)).c_str(), &vec->Data[0].x, &vec->Data[0].y,
                       vec->Data.size(), 0, vec->Offset, 2 * sizeof(float));

    ImPlot::EndPlot();
  }
}

inline void HelpMarker(const char *desc) {
  ImGui::TextDisabled("(?)");
  if (ImGui::IsItemHovered()) {
    ImGui::BeginTooltip();
    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
    ImGui::TextUnformatted(desc);
    ImGui::PopTextWrapPos();
    ImGui::EndTooltip();
  }
}

// inline void displayVector3d(const char *vecName, Eigen::Vector3d vec) {
//   ImGui::TableNextRow();
//   ImGui::TableSetColumnIndex(0);
//   ImGui::Text("%s", vecName);
//   for (int i = 0; i < 3; ++i) {
//     ImGui::TableSetColumnIndex(i + 1);
//     ImGui::Text("%f", vec(i));
//   }
// }

inline void displayArray3d(const char *arrName, double arr[3], const char *helpText = "") {
  ImGui::TableNextRow();
  ImGui::TableSetColumnIndex(0);
  ImGui::Text("%s", arrName);
  ImGui::SameLine();
  HelpMarker(helpText);
  for (int i = 0; i < 3; ++i) {
    ImGui::TableSetColumnIndex(i + 1);
    ImGui::Text("%f", arr[i]);
  }
}

inline void displayArrayNd(const char *arrName, const std::vector<double> &arr,
                           const char *helpText = "") {
  ImGui::TableNextRow();
  ImGui::TableSetColumnIndex(0);
  ImGui::Text("%s", arrName);
  ImGui::SameLine();
  HelpMarker(helpText);
  for (size_t i = 0; i < arr.size(); ++i) {
    ImGui::TableSetColumnIndex(i + 1);
    ImGui::Text("%f", arr[i]);
  }
}

inline void displayArray3b(const char *arrName, bool arr[3], const char *helpText = "") {
  ImGui::TableNextRow();
  ImGui::TableSetColumnIndex(0);
  ImGui::Text("%s", arrName);
  ImGui::SameLine();
  HelpMarker(helpText);
  for (int i = 0; i < NUM_CHANS; ++i) {
    ImGui::TableSetColumnIndex(i + 1);
    ImGui::Text("%d", arr[i]);
  }
}
