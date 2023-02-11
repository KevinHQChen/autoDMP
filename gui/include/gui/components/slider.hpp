#pragma once

#include "imgui.h"
#include <functional>
#include <string>
#include <type_traits>

namespace gui {

template <typename T> class Slider {
  std::string label_;
  T min_;
  T max_;
  std::string format_;
  int numSliders_;
  ImVec2 sliderSize_;
  T *values_;
  bool horizontal_;
  std::function<void()> callback_;

public:
  Slider(std::string label, T min, T max, std::string format = "%d", int numSliders = 1,
         ImVec2 sliderSize = ImVec2(36, 200), bool horizontal = true,
         std::function<void()> callback = nullptr)
      : label_(label), min_(min), max_(max), format_(format), numSliders_(numSliders),
        sliderSize_(sliderSize), values_(new T[numSliders]), horizontal_(horizontal),
        callback_(callback) {
    for (int i = 0; i < numSliders_; ++i)
      values_[i] = min;
  }

  template <typename U = T>
  typename std::enable_if<std::is_integral<U>::value, void> renderIntSlider(std::string label,
                                                                            int idx) {
    if (horizontal_)
      ImGui::SliderInt(label.c_str(), (int *)&values_[idx], min_, max_, format_.c_str());
    else
      ImGui::VSliderInt(label.c_str(), sliderSize_, (int *)&values_[idx], min_, max_,
                        format_.c_str());
  }

  template <typename U = T>
  typename std::enable_if<std::is_floating_point<U>::value, void>
  renderFloatSlider(std::string label, int idx) {
    if (horizontal_)
      ImGui::SliderFloat(label.c_str(), (float *)&values_[idx], min_, max_, format_.c_str());
    else
      ImGui::VSliderFloat(label.c_str(), sliderSize_, (float *)&values_[idx], min_, max_,
                          format_.c_str());
  }

  void render() {
    ImGui::BeginGroup();
    ImGui::Text("%s", label_.c_str());
    for (int i = 0; i < numSliders_; ++i) {
      if (i != 0 && !horizontal_)
        ImGui::SameLine();
      std::string label = "##" + label_ + std::to_string(i);
      if (std::is_integral<T>::value)
        renderIntSlider<T>(label, i);
      else if (std::is_floating_point<T>::value)
        renderFloatSlider<T>(label, i);
    }
    ImGui::EndGroup();
    if (callback_)
      callback_();
  }

  T get(int index) const { return values_[index]; }

  void set(int index, T value) { values_[index] = value; }

  T *get() const { return values_; }
};

} // namespace gui
