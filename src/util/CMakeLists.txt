# Utility library for autoDMP
set(HEADER_LIST "../../include/util/util.hpp")

# Make an automatic library - will be static or dynamic based on user setting
add_library(util util.cpp ${HEADER_LIST})   # target is created

find_package(fmt CONFIG)
find_package(spdlog CONFIG)
find_package(CLI11 CONFIG REQUIRED)
find_package(toml11 CONFIG REQUIRED)
find_package(tsl-ordered-map CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)
find_package(SDL2 CONFIG REQUIRED)

find_package(OpenCV REQUIRED)
if (${OpenCV_FOUND})
    message("OpenCV include dirs: ${OpenCV_INCLUDE_DIRS}")
    include_directories(${OpenCV_INCLUDE_DIRS})
    message("OpenCV link libraries: ${OpenCV_LIBS}")
else()
    message("Couldn't find OpenCV")
endif()

find_package(OpenGL REQUIRED)
if (OPENGL_FOUND)
    message("OpenGL include dir: ${OPENGL_INCLUDE_DIR}")
    include_directories( ${OPENGL_INCLUDE_DIR} )
    message("OpenGL link libraries: ${OPENGL_opengl_LIBRARY}")
else (OPENGL_FOUND)
    message("Couldn't find OpenGL")
endif()

find_package(Boost REQUIRED COMPONENTS program_options)
if (${Boost_FOUND})
    include_directories( ${Boost_INCLUDE_DIR} )
else()
    message("Couldn't find Boost")
endif()

target_include_directories(util
  PRIVATE "${CMAKE_BINARY_DIR}/configured_files/include"
)

target_link_libraries(
  util
  PUBLIC project_options project_warnings
  PRIVATE CLI11::CLI11
  PRIVATE fmt::fmt fmt::fmt-header-only
  PRIVATE spdlog::spdlog spdlog::spdlog_header_only
  PRIVATE SDL2::SDL2 SDL2::SDL2main SDL2::SDL2-static
  PRIVATE toml11::toml11
  PRIVATE tsl::ordered_map
  PRIVATE imgui::imgui
  PUBLIC ${OpenCV_LIBS}
  PUBLIC ${OPENGL_opengl_LIBRARY}
  PUBLIC Boost::program_options
)

# We need this directory, and users of our library will need it too
target_include_directories(util
    PUBLIC ../../include
    PUBLIC ${OpenCV_INCLUDE_DIRS}
    PUBLIC ${Boost_INCLUDE_DIR}
    PUBLIC ${OPENGL_INCLUDE_DIR}
)
