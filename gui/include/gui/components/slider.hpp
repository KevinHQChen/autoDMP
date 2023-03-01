#pragma once

#include "imgui.h"
#include <cassert>
#include <functional>
#include <string>
#include <type_traits>

namespace gui {

template <typename T> class Slider {
  std::string label_;
  T min_, max_;
  T *value_{nullptr};
  std::string format_;
  ImVec2 sliderSize_;
  bool horizontal_;
  std::function<void()> callback_;

public:
  Slider(std::string label, T min, T max, T *value, std::string format = "%d",
         ImVec2 sliderSize = ImVec2(36, 200), bool horizontal = true,
         std::function<void()> callback = nullptr)
      : label_(label), min_(min), max_(max), format_(format), value_(value),
        sliderSize_(sliderSize), horizontal_(horizontal), callback_(callback) {
    if (value_ == nullptr)
      value_ = new T(min_);
    assert(value_ != nullptr);
  }

  template <typename U = T>
  typename std::enable_if<std::is_integral<U>::value, void> renderIntSlider(std::string label) {
    if (horizontal_)
      ImGui::SliderInt(label.c_str(), (int *)value_, min_, max_, format_.c_str());
    else
      ImGui::VSliderInt(label.c_str(), sliderSize_, (int *)value_, min_, max_, format_.c_str());
  }

  template <typename U = T>
  typename std::enable_if<std::is_floating_point<U>::value, void>
  renderFloatSlider(std::string label) {
    if (horizontal_)
      ImGui::SliderFloat(label.c_str(), (float *)value_, min_, max_);
    else
      ImGui::VSliderFloat(label.c_str(), sliderSize_, (float *)value_, min_, max_);
  }

  void render() {
    ImGui::BeginGroup();
    // TODO set slider size:
    // ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
    if (std::is_integral<T>::value)
      renderIntSlider<T>(label_);
    else if (std::is_floating_point<T>::value)
      renderFloatSlider<T>(label_);
    ImGui::EndGroup();
    if (callback_)
      callback_();
  }

  T get() const { return *value_; }

  void set(T value) { *value_ = value; }
};

template <typename T> class SliderArray {
  std::string label_;
  T min_, max_;
  int numSliders_;
  std::vector<T> *values_{nullptr};
  std::string format_;
  ImVec2 sliderSize_;
  bool horizontal_;
  std::function<void()> callback_;

public:
  SliderArray(std::string label, T min, T max, std::vector<T> *values = nullptr, int numSliders = 1,
              std::string format = "%d", ImVec2 sliderSize = ImVec2(36, 200),
              bool horizontal = true, std::function<void()> callback = nullptr)
      : label_(label), min_(min), max_(max), values_(values), format_(format),
        sliderSize_(sliderSize), horizontal_(horizontal), callback_(callback) {
    if (values_ == nullptr) {
      values_ = new std::vector<T>(numSliders, min_);
      numSliders_ = numSliders;
    } else
      numSliders_ = values_->size();
    assert(values_->size() == numSliders_);
  }

  template <typename U = T>
  typename std::enable_if<std::is_integral<U>::value, void> renderIntSlider(std::string label,
                                                                            int idx) {
    if (horizontal_)
      ImGui::SliderInt(label.c_str(), (int *)&(*values_)[idx], min_, max_, format_.c_str());
    else
      ImGui::VSliderInt(label.c_str(), sliderSize_, (int *)&(*values_)[idx], min_, max_,
                        format_.c_str());
  }

  template <typename U = T>
  typename std::enable_if<std::is_floating_point<U>::value, void>
  renderFloatSlider(std::string label, int idx) {
    if (horizontal_)
      ImGui::SliderFloat(label.c_str(), (float *)&(*values_)[idx], min_, max_);
    else
      ImGui::VSliderFloat(label.c_str(), sliderSize_, (float *)&(*values_)[idx], min_, max_);
  }

  void render() {
    ImGui::BeginGroup();
    ImGui::Text("%s", label_.c_str());
    for (int i = 0; i < numSliders_; ++i) {
      std::string label = "##" + label_ + std::to_string(i);
      if (std::is_integral<T>::value)
        renderIntSlider<T>(label, i);
      else if (std::is_floating_point<T>::value)
        renderFloatSlider<T>(label, i);
      if (!horizontal_)
        ImGui::SameLine();
    }
    ImGui::EndGroup();
    if (callback_)
      callback_();
  }

  T get(int index) const { return (*values_)[index]; }

  std::vector<T> get() const { return *values_; }

  void set(T value, int index = 0) { (*values_)[index] = value; }
};

} // namespace gui
