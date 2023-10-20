cmake_minimum_required(VERSION 3.23)

# ---- Project ----

# Setup vcpkg before defining the project
set(CMAKE_TOOLCHAIN_FILE ~/vcpkg/scripts/buildsystems/vcpkg.cmake CACHE FILEPATH "Path to toolchain")

# Set the project name and language
project(
  AutoDMP
  VERSION 0.1.0
  DESCRIPTION "An automatic droplet microfluidics platform"
  HOMEPAGE_URL "https://github.com/KevinHQChen/autoDMP"
  LANGUAGES CXX C)

# ---- Features ----

include(cmake/features.cmake)

# ---- Source files ----

add_subdirectory(./apps)
add_subdirectory(./util)
add_subdirectory(./cam)
add_subdirectory(./imcap)
add_subdirectory(./improc)
add_subdirectory(./ctrl)
add_subdirectory(./gui)
add_subdirectory(./pump)
add_subdirectory(./external)
add_subdirectory(./pump_sinewave)
