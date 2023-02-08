#pragma once

#include "imgui.h"
#include <functional>
#include <string>

namespace gui {

class Button {
  std::vector<std::string> labels_;
  int numButtons_;
  ImVec2 buttonSize_;
  std::function<void()> callback_;

public:
  Button(std::vector<std::string> labels, std::function<void()> callback = nullptr,
         int numButtons = 1, ImVec2 buttonSize = ImVec2(0, 0))
      : labels_(labels), numButtons_(numButtons), buttonSize_(buttonSize), callback_(callback) {}

  void render() {
    ImGui::BeginGroup();
    for (int i = 0; i < numButtons_; i++) {
      if (i != 0)
        ImGui::SameLine();
      bool pressed = ImGui::Button(labels_[i].c_str(), buttonSize_);
      if (pressed && callback_)
        callback_();
    }
    ImGui::EndGroup();
  }

  void setLabel(int idx, std::string label) { labels_[idx] = label; }
};

} // namespace gui
