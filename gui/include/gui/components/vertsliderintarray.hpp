#pragma once

#include "imgui.h"
#include <functional>
#include <string>

namespace gui {

class VertSliderIntArray {
  std::string label_;
  int size_;
  ImVec2 sliderSize_;
  std::string format_;
  int *values_;
  int min_;
  int max_;
  std::function<void()> callback_;

public:
  VertSliderIntArray(std::string label, int size, int min, int max,
                     ImVec2 sliderSize = ImVec2(36, 200), std::string format = "%d", std::function<void()> callback = nullptr)
      : label_(label), size_(size), sliderSize_(sliderSize), format_(format),
        values_(new int[size]), min_(min), max_(max), callback_(callback) {
    for (int i = 0; i < size; ++i) {
      values_[i] = 0;
    }
  }

  void render() {
    ImGui::BeginGroup();
    for (int i = 0; i < size_; ++i) {
      std::string label = "##" + label_ + std::to_string(i);
      ImGui::VSliderInt(label.c_str(), sliderSize_, &values_[i], min_, max_, format_.c_str());
      ImGui::SameLine();
    }
    ImGui::Text("%s", label_.c_str());
    ImGui::EndGroup();
    callback_();
  }

  int *getValues() { return values_; }
};

} // namespace gui
