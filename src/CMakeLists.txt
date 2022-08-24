# add_subdirectory(sdl)
add_subdirectory(fltk)

set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0)
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(fmt CONFIG)
find_package(spdlog CONFIG)
find_package(CLI11 CONFIG REQUIRED)

# Generic test that uses conan libs
add_executable(intro main.cpp)
target_link_libraries(
  intro
  PUBLIC project_options project_warnings
  PRIVATE CLI11::CLI11 fmt::fmt fmt::fmt-header-only spdlog::spdlog spdlog::spdlog_header_only)

target_include_directories(intro PRIVATE "${CMAKE_BINARY_DIR}/configured_files/include")
