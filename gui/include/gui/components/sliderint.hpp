#pragma once

#include "imgui.h"
#include <functional>
#include <string>

namespace gui {

class SliderInt {
  std::string label_;
  int min_;
  int max_;
  int value_;
  std::function<void()> callback_;

public:
  SliderInt(std::string label, int min, int max)
      : label_(label), min_(min), max_(max), value_(min) {}

  void render() {
    ImGui::SliderInt(label_.c_str(), &value_, min_, max_);
    ImGui::Text("%s", label_.c_str());
  }

  int get() const { return value_; }
};

} // namespace gui
