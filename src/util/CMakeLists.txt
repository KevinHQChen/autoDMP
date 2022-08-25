# Utility library for autoDMP
# Note that headers are optional, and do not affect add_library, but they will
# not show up in IDEs unless they are listed in add_library.

set(HEADER_LIST "../../include/util/util.hpp")

# Make an automatic library - will be static or dynamic based on user setting
add_library(util util.cpp ${HEADER_LIST})   # target is created

# one of these libraries needs pthread
set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0)
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(fmt CONFIG)
find_package(spdlog CONFIG)
find_package(CLI11 CONFIG REQUIRED)

find_package(OpenCV REQUIRED)
if (${OpenCV_FOUND})
    include_directories(${OpenCV_INCLUDE_DIRS})
else()
    message("Couldn't find OpenCV")
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
  PRIVATE CLI11::CLI11 fmt::fmt fmt::fmt-header-only spdlog::spdlog spdlog::spdlog_header_only
  PUBLIC ${OpenCV_LIBS}
  PUBLIC Boost::program_options
)

# We need this directory, and users of our library will need it too
target_include_directories(util
    PUBLIC ../../include
    PUBLIC ../../external/eigen
    PUBLIC ${OpenCV_INCLUDE_DIRS}
    PUBLIC ${Boost_INCLUDE_DIR}
    )
