#pragma once

#include "imgui.h"
#include <cassert>
#include <functional>
#include <iostream>
#include <string>
#include <type_traits>

namespace gui {

template <typename T> class Slider {
  std::string label_;
  T min_, max_;
  T *value_;
  bool internalValue_{false};
  std::string format_;
  ImVec2 sliderSize_;
  bool horizontal_;
  std::function<void()> callback_;

public:
  Slider(
      std::string label, T min, T max, T *value = nullptr, std::string format = "%d",
      ImVec2 sliderSize = ImVec2(36, 200), bool horizontal = true,
      std::function<void()> callback = []() {})
      : label_(label), min_(min), max_(max), format_(format), value_(value),
        sliderSize_(sliderSize), horizontal_(horizontal), callback_(callback) {
    if (value_ == nullptr) {
      value_ = new T(min_);
      internalValue_ = true;
    }
  }

  ~Slider() {
    if (internalValue_)
      delete value_;
  }

  void render() {
    ImGui::BeginGroup();
    // TODO set slider size:
    // ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
    if (horizontal_) {
      if constexpr (std::is_same<T, int>::value)
        ImGui::SliderInt(label_.c_str(), value_, min_, max_, format_.c_str());
      else if constexpr (std::is_same<T, float>::value)
        ImGui::SliderFloat(label_.c_str(), value_, min_, max_);
    } else {
      if constexpr (std::is_same<T, int>::value)
        ImGui::VSliderInt(label_.c_str(), sliderSize_, value_, min_, max_, format_.c_str());
      else if constexpr (std::is_same<T, float>::value)
        ImGui::VSliderFloat(label_.c_str(), sliderSize_, value_, min_, max_);
    }
    ImGui::EndGroup();
    callback_();
  }

  T get() const { return *value_; }

  void set(T value) { *value_ = value; }
};

template <typename T> class SliderArray {
  std::string label_;
  T min_, max_;
  int numSliders_;
  std::vector<T> *values_;
  bool internalValues_{false};
  std::string format_;
  ImVec2 sliderSize_;
  bool horizontal_;
  std::function<void()> callback_;

public:
  SliderArray(
      std::string label, T min, T max, std::vector<T> *values = nullptr, int numSliders = 1,
      std::string format = "%d", ImVec2 sliderSize = ImVec2(36, 200), bool horizontal = true,
      std::function<void()> callback = []() {})
      : label_(label), min_(min), max_(max), values_(values), format_(format),
        sliderSize_(sliderSize), horizontal_(horizontal), callback_(callback) {
    if (values_ == nullptr) {
      values_ = new std::vector<T>(numSliders, min_);
      internalValues_ = true;
      numSliders_ = numSliders;
    } else
      numSliders_ = values_->size();
    assert(values_->size() == numSliders_);
  }

  ~SliderArray() {
    if (internalValues_)
      delete values_;
  }

  void render() {
    ImGui::BeginGroup();
    ImGui::Text("%s", label_.c_str());
    for (int i = 0; i < numSliders_; ++i) {
      std::string label = "##" + label_ + std::to_string(i);
      if (horizontal_) {
        if constexpr (std::is_same<T, int>::value)
          ImGui::SliderInt(label.c_str(), &(*values_)[i], min_, max_, format_.c_str());
        else if constexpr (std::is_same<T, float>::value)
          ImGui::SliderFloat(label.c_str(), &(*values_)[i], min_, max_);
      } else {
        if constexpr (std::is_same<T, int>::value)
          ImGui::VSliderInt(label.c_str(), sliderSize_, &(*values_)[i], min_, max_,
                            format_.c_str());
        else if constexpr (std::is_same<T, float>::value)
          ImGui::VSliderFloat(label.c_str(), sliderSize_, &(*values_)[i], min_, max_);
        ImGui::SameLine();
      }
    }
    ImGui::EndGroup();
    callback_();
  }

  T get(int index) const { return (*values_)[index]; }

  std::vector<T> get() const { return *values_; }

  void set(T value, int index = 0) { (*values_)[index] = value; }

  void setMax(int max) {
    if (max != max_) {
      max_ = max;
      for (auto &v : *values_)
        v = std::min(v, max_);
    }
  }
};

} // namespace gui
