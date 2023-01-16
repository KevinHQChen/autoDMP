#pragma once

#include "imgui.h"
#include <functional>
#include <string>

namespace gui {

class CheckboxArray {
  std::string label_;
  int size_;
  bool* enabled_;

  std::function<void()> callback_;
public:
  CheckboxArray(std::string label, int size)
      : label_(label), size_(size), enabled_(new bool[size_]) {
  for (int i = 0; i < size; ++i)
    enabled_[i] = false;
  }

  CheckboxArray(std::string label, int size, std::function<void()> callback)
      : label_(label), size_(size), enabled_(new bool[size_]), callback_(callback) {
  for (int i = 0; i < size; ++i)
    enabled_[i] = false;
  }

  void render() {
    ImGui::BeginGroup();

    for (int i = 0; i < size_; ++i) {
      ImGui::Checkbox((label_ + std::to_string(i)).c_str(), &enabled_[i]);
      ImGui::SameLine();
    }
    ImGui::Text("");

    ImGui::EndGroup();
  }

  bool* get() {
    return enabled_;
  }
};

} // namespace gui
