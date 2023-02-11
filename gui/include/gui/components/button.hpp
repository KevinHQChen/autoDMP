#pragma once

#include "imgui.h"
#include <functional>
#include <string>

namespace gui {

class Button {
  std::vector<std::string> labels_;
  int numButtons_;
  ImVec2 buttonSize_;
  bool *values_;
  std::function<void()> callback_;

public:
  Button(std::string label, std::function<void()> callback = nullptr,
         ImVec2 buttonSize = ImVec2(0, 0))
      : numButtons_(1), buttonSize_(buttonSize), callback_(callback) {
    labels_.push_back(label);
    values_ = new bool[numButtons_];
  }

  Button(std::vector<std::string> labels, std::function<void()> callback = nullptr,
         int numButtons = 1, ImVec2 buttonSize = ImVec2(0, 0))
      : labels_(labels), numButtons_(numButtons), buttonSize_(buttonSize),
        values_(new bool[numButtons_]), callback_(callback) {}

  void render() {
    ImGui::BeginGroup();
    for (int i = 0; i < numButtons_; i++) {
      if (i != 0)
        ImGui::SameLine();
      values_[i] = ImGui::Button(labels_[i].c_str(), buttonSize_);
    }
    ImGui::EndGroup();
    if (callback_)
      callback_();
  }

  void setLabel(int idx, std::string label) { labels_[idx] = label; }

  bool get(int index) const { return values_[index]; }

  bool *get() const { return values_; }
};

} // namespace gui
