#pragma once

#include <functional>
#include <memory>
#include "util/util.hpp"

// include/gui/windows/window.hpp
// base class for all windows
// defines a pure virtual render method
// defines a registerCallback method to register callback functions for child components
// defines a handleEvent method called in the render method to handle all registered events for
// child components at once
namespace gui {
class Window {
public:
  Window() = default;
  virtual ~Window() = default;

  virtual void render() = 0;

protected:
  void registerCallback(std::function<void()> callback);

  void handleEvent();

private:
  std::vector<std::function<void()>> callbacks_;
};
} // namespace gui
