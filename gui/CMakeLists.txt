add_library(gui
  "./src/gui/gui.cpp"
  "./src/gui/windows/window.cpp"
  "./src/gui/windows/improc_window.cpp"
  "./src/gui/windows/pump_window.cpp"
  "./src/gui/windows/ctrl_window.cpp"
  "./src/gui/windows/plot_window.cpp")

# includes
include(FetchContent)
Set(FETCHCONTENT_QUIET FALSE)
FetchContent_Declare(
    imgui_bundle
    GIT_REPOSITORY https://github.com/pthom/imgui_bundle.git
    GIT_PROGRESS TRUE
    # Enter the desired git tag below
    GIT_TAG "v0.8.5"
)
FetchContent_MakeAvailable(imgui_bundle)

set(gui_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  CACHE STRING "")
target_include_directories(gui PUBLIC
  "$<BUILD_INTERFACE:${gui_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:./${CMAKE_INSTALL_INCLUDEDIR}>")

# link dependencies
target_link_libraries(gui PRIVATE
  util
  ctrl
  supervisor
  improc
  imcap
  cam
  pump
  imgui_bundle
)

# tests
if(BUILD_TESTS)
  add_subdirectory("./test")
endif()
