#pragma once

#include "imgui.h"
#include <functional>
#include <string>

namespace gui {

class SliderFloatArray {
public:
  SliderFloatArray(std::string label, int size, float min, float max,
                   std::function<void()> callback)
      : label_(label), size_(size), values_(new float[size]), min_(min), max_(max),
        callback_(callback) {
    for (int i = 0; i < size; ++i) {
      values_[i] = 0.0f;
    }
  }

  void render() {
    ImGui::BeginGroup();
    for (int i = 0; i < size_; ++i) {
      ImGui::SetNextItemWidth(ImGui::GetWindowWidth() / size_ / 2);
      std::string label = "##" + label_ + std::to_string(i);
      ImGui::SliderFloat(label.c_str(), &values_[i], min_, max_);
      ImGui::SameLine();
    }
    ImGui::Text("%s", label_.c_str());
    ImGui::EndGroup();
  }

private:
  std::string label_;
  int size_;
  float *values_;
  float min_;
  float max_;
  std::function<void()> callback_;
};

} // namespace gui
