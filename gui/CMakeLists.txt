add_library(gui
  "./src/gui/gui.cpp"
  "./src/gui/windows/window.cpp"
  # "./src/gui/windows/imcap_window.cpp"
  # "./src/gui/windows/improc_setup_window.cpp"
  # "./src/gui/windows/improc_window.cpp"
  # "./src/gui/windows/pump_window.cpp"
  "./src/gui/windows/sysid_setup_window.cpp")
  # "./src/gui/windows/sysid_window.cpp"
  # "./src/gui/windows/ctrl_setup_window.cpp"
  # "./src/gui/windows/ctrl_window.cpp")

target_link_libraries(gui PRIVATE project_options project_warnings)

# includes
set(gui_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  CACHE STRING "")
target_include_directories(gui PUBLIC
  "$<BUILD_INTERFACE:${gui_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:./${CMAKE_INSTALL_INCLUDEDIR}>")

# find dependencies
set(gui_DEPENDENCIES_CONFIGURED
    imgui
    implot
    OpenGL
    GLEW
    # glfw3
    # SDL2
    CACHE STRING "")
foreach (dep ${gui_DEPENDENCIES_CONFIGURED})
  find_package(${dep} REQUIRED)
  if (DEFINED ${dep}_INCLUDE_DIRS)
    target_include_directories(gui PUBLIC ${${dep}_INCLUDE_DIRS})
  endif ()
endforeach ()

# link dependencies
target_link_system_libraries(gui PRIVATE
  util
  ctrl
  improc
  imcap
  cam
  pump
  imgui::imgui
  implot::implot
  ${OPENGL_opengl_LIBRARY}
  GLEW::GLEW
  glfw
  # SDL2::SDL2 SDL2::SDL2main SDL2::SDL2-static
)

# tests
if(FEATURE_TESTS)
  add_subdirectory("./test")
endif()
