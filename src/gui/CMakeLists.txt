# GUI library for autoDMP
set(HEADER_LIST "../../include/gui/gui.hpp")

# Make an automatic library - will be static or dynamic based on user setting
add_library(gui gui.cpp ${HEADER_LIST})   # target is created

find_package(imgui CONFIG REQUIRED)
find_package(SDL2 CONFIG REQUIRED)

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
    message("GLEW link libraries: ${GLEW_LIBRARIES}")
else ()
    message("Couldn't find GLEW")
endif()

find_package(glfw3 CONFIG REQUIRED)
# if (GLFW_FOUND)
#     message("GLFW3 include dir: ${GLFW3_INCLUDE_DIR}")
#     include_directories( ${GLFW3_INCLUDE_DIR} )
#     message("GLFW3 link libraries: ${GLFW3_LIBRARIES}")
# else ()
#     message("Couldn't find GLFW3")
# endif()

target_link_libraries(
  gui
  PUBLIC util
  PRIVATE imgui::imgui
  PRIVATE SDL2::SDL2 SDL2::SDL2main SDL2::SDL2-static
  PRIVATE ${OPENGL_opengl_LIBRARY}
  # PRIVATE ${OPENGL_glu_LIBRARY}
  PRIVATE ${GLEW_LIBRARIES}
  PRIVATE glfw
)

# We need this directory, and users of our library will need it too
target_include_directories(util
    PUBLIC ../../include
    PUBLIC ${OPENGL_INCLUDE_DIR}
)
