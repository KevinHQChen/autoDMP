# GUI library for autoDMP
set(HEADER_LIST "../../include/gui/gui.hpp"
    "../../include/gui/guiframe.hpp")

# Make an automatic library - will be static or dynamic based on user setting
add_library(gui gui.cpp ${HEADER_LIST})   # target is created

# include(FetchContent)
# FetchContent_Declare(_imguiwrap URL https://github.com/kfsone/imguiwrap/archive/refs/tags/v1.2.2.tar.gz)
# FetchContent_Declare(_imguiwrap GIT_REPOSITORY https://github.com/KevinHQChen/imguiwrap.git
#   GIT_TAG origin/main)
# FetchContent_MakeAvailable(_imguiwrap)
# message("imguiwrap vars: ${imguiwrap_SOURCE_DIR}")
# include_directories(${imguiwrap_SOURCE_DIR})

find_package(imgui CONFIG REQUIRED)
# find_package(SDL2 CONFIG REQUIRED)

find_package(OpenGL REQUIRED)
if (OPENGL_FOUND)
    message("OpenGL include dir: ${OPENGL_INCLUDE_DIR}")
    include_directories( ${OPENGL_INCLUDE_DIR} )
    message("OpenGL link libraries: ${OPENGL_opengl_LIBRARY}")
else ()
    message("Couldn't find OpenGL")
endif()

find_package(GLEW REQUIRED)
if (GLEW_FOUND)
    message("GLEW include dir: ${GLEW_INCLUDE_DIRS}")
    include_directories( ${GLEW_INCLUDE_DIRS} )
else ()
    message("Couldn't find GLEW")
endif()

# find_package(glfw3 CONFIG REQUIRED)
# if (GLFW3_FOUND)
#     message("GLFW3 include dir: ${GLFW3_INCLUDE_DIR}")
#     include_directories( ${GLFW3_INCLUDE_DIR} )
#     message("GLFW3 link libraries: ${GLFW3_LIBRARIES}")
# else ()
#     message("Couldn't find GLFW3")
# endif()

target_link_libraries(
  gui
  PUBLIC util
  PUBLIC imcap
  # PUBLIC imguiwrap
  PUBLIC imgui::imgui
  # PRIVATE SDL2::SDL2 SDL2::SDL2main SDL2::SDL2-static
  PRIVATE ${OPENGL_opengl_LIBRARY}
  # PRIVATE ${OPENGL_glu_LIBRARY}
  PRIVATE GLEW::GLEW
  PRIVATE glfw
)

# We need this directory, and users of our library will need it too
target_include_directories(util
    PUBLIC ../../include
    PUBLIC ${OPENGL_INCLUDE_DIR}
    PUBLIC ${imguiwrap_SOURCE_DIR}
)
