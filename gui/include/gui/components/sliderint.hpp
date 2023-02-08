#pragma once

#include "imgui.h"
#include <functional>
#include <string>

namespace gui {

class SliderInt {
  std::string label_;
  int min_;
  int max_;
  std::string format_;
  int numSliders_;
  ImVec2 sliderSize_;
  int *values_;
  bool horizontal_;
  std::function<void()> callback_;

public:
  SliderInt(std::string label, int min, int max, std::string format = "%d", int numSliders = 1,
            ImVec2 sliderSize = ImVec2(36, 200), bool horizontal = true,
            std::function<void()> callback = nullptr)
      : label_(label), min_(min), max_(max), format_(format), numSliders_(numSliders),
        sliderSize_(sliderSize), values_(new int[numSliders]), horizontal_(horizontal),
        callback_(callback) {
    for (int i = 0; i < numSliders_; ++i)
      values_[i] = min;
  }

  void render() {
    ImGui::BeginGroup();
    ImGui::Text("%s", label_.c_str());
    for (int i = 0; i < numSliders_; ++i) {
      if (i != 0 && !horizontal_)
        ImGui::SameLine();
      std::string label = "##" + label_ + std::to_string(i);
      if (horizontal_)
        ImGui::SliderInt(label.c_str(), &values_[i], min_, max_, format_.c_str());
      else
        ImGui::VSliderInt(label.c_str(), sliderSize_, &values_[i], min_, max_, format_.c_str());
    }
    ImGui::EndGroup();
    if (callback_)
      callback_();
  }

  int get(int index) const { return values_[index]; }

  int *get() const { return values_; }
};

} // namespace gui
