cmake_minimum_required(VERSION 3.21)

# ensure tooling and projects use the same version of C++
set(CMAKE_CXX_STANDARD 20) # C++20
set(CMAKE_CXX_STANDARD_REQUIRED ON) # don't allow compiler to fall back to older standard
set(CMAKE_CXX_EXTENSIONS OFF) # disable compiler-specific extensions (force strict adherence to standard)

# ---- Project ----

# Add project_options (update version number regularly - manually for now)
include(FetchContent)
FetchContent_Declare(_project_options URL https://github.com/aminya/project_options/archive/refs/tags/v0.26.3.zip)
FetchContent_MakeAvailable(_project_options)
include(${_project_options_SOURCE_DIR}/Index.cmake)

# Setup vcpkg before defining the project
run_vcpkg()

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
# add_subdirectory(./gui)

# ---- Package ----

if(CMAKE_SKIP_INSTALL_RULES)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  return()
endif()

# TODO where is this set?
set(GIT_SHA
    "Unknown"
    CACHE STRING "SHA this build was generated from")
string(
  SUBSTRING "${GIT_SHA}"
            0
            8
            GIT_SHORT_SHA)

# the variables set using CACHE STRING "" are passed to package_project

package_project(
  TARGETS
  autoDMP
  # util
  # cam
  # gui
  # imcap
  # improc
  project_options
  project_warnings
  INTERFACE_DEPENDENCIES_CONFIGURED
  # ${my_header_lib_DEPENDENCIES_CONFIGURED}
  INTERFACE_INCLUDES
  # ${my_header_lib_INCLUDE_DIR}
  PUBLIC_INCLUDES
  # ${util_INCLUDE_DIR}
  # ${cam_INCLUDE_DIR}
  # ${gui_INCLUDE_DIR}
  # ${imcap_INCLUDE_DIR}
  # ${improc_INCLUDE_DIR}
)

set(CPACK_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}-${GIT_SHORT_SHA}-${CMAKE_SYSTEM_NAME}-${CMAKE_BUILD_TYPE}-${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}"
)

include(CPack)
