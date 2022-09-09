cmake_minimum_required(VERSION 3.21)

# Not ideal to use this global variable, but necessary to make sure
# that tooling and projects use the same version
set(CMAKE_CXX_STANDARD 20)

# strongly encouraged to enable this globally to avoid conflicts between
# -Wpedantic being enabled and -std=c++20 and -std=gnu++20 for example
# when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)

# Note: by default ENABLE_DEVELOPER_MODE is True
# This means that all analysis (sanitizers, static analysis)
# is enabled and all warnings are treated as errors
# if you want to switch this behavior, change TRUE to FALSE
set(ENABLE_DEVELOPER_MODE
    FALSE
    CACHE BOOL "Enable 'developer mode'")
set(OPT_WARNINGS_AS_ERRORS_DEVELOPER_DEFAULT FALSE) # disable warnings_as_errors in developer mode (for now)

# Add project_options (https://github.com/cpp-best-practices/project_options)
include(FetchContent)
FetchContent_Declare(_project_options URL https://github.com/aminya/project_options/archive/refs/tags/v0.24.1.zip)
FetchContent_MakeAvailable(_project_options)
include(${_project_options_SOURCE_DIR}/Index.cmake)

# Setup vcpkg - should be called before defining project()
run_vcpkg()

# Set the project name and language
project(
  autoDMP
  VERSION 0.1.0
  DESCRIPTION ""
  HOMEPAGE_URL "https://github.com/KevinHQChen/autoDMP"
  LANGUAGES CXX C)

if(PROJECT_IS_TOP_LEVEL)  # this variable is set in project()
  # Consider the CTest module, which creates targets and options!
  # Only needed if you want to enable submissions to a CDash server.
  include(CTest)
endif()

set(GIT_SHA
    "Unknown"
    CACHE STRING "SHA this build was generated from")
string(
  SUBSTRING "${GIT_SHA}"
            0
            8
            GIT_SHORT_SHA)

if(GENERATOR_IS_MULTI_CONFIG AND NOT CMAKE_BUILD_TYPE)
  # Make sure that all supported configuration types have their
  # associated conan packages available. You can reduce this
  # list to only the configuration types you use, but only if one
  # is not forced-set on the command line for VS
  message(TRACE "Setting up multi-config build types")
  set(CMAKE_CONFIGURATION_TYPES
      Debug
      Release
      RelWithDebInfo
      MinSizeRel
      CACHE STRING "Enabled build types" FORCE)
else()
  message(TRACE "User chose a specific build type, so we are using that")
  set(CMAKE_CONFIGURATION_TYPES
      ${CMAKE_BUILD_TYPE}
      CACHE STRING "Enabled build types" FORCE)
endif()

# Add dynamic cmake options
include(${_project_options_SOURCE_DIR}/src/DynamicProjectOptions.cmake)

# wrapper macro for project_options() that sets recommended defaults,
# provides user/developer modes and full GUI support for choosing options at configure time
#
# default options:
# ENABLE_DEVELOPER_MODE: ON
#  * WARNINGS_AS_ERRORS: ON
#  * ENABLE_SANITIZER_ADDRESS: ON
#  * ENABLE_CLANG_TIDY: ON for Ninja/Makefiles
#  * ENABLE_SANITIZER_UNDEFINED: ON for Compilers that support it
#  * ENABLE_CPPCHECK: ON for Ninja/Makefiles
#
# Any default can be overridden
# set(<feature_name>_DEFAULT <value>) - set default for both user and developer modes
# set(<feature_name>_DEVELOPER_DEFAULT <value>) - set default for developer mode
# set(<feature_name>_USER_DEFAULT <value>) - set default for user mode
dynamic_project_options(
  # Note: PCH is disabled by default in developer mode because these headers become
  # globally included and they can mask other errors
  # common headers to pre-compile
  PCH_HEADERS
  <vector>
  <string>
  # MSVC_WARNINGS    # Override the defaults for the MSVC warnings
  # CLANG_WARNINGS   # Override the defaults for the CLANG warnings
  # GCC_WARNINGS     # Override the defaults for the GCC warnings
  CPPCHECK_OPTIONS
  --enable=style,performance,warning,portability
  --inline-suppr
  # We cannot act on a bug/missing feature of cppcheck
  --suppress=cppcheckError
  --suppress=internalAstError
  # if a file does not have an internalAstError, we get an unmatchedSuppression error
  --suppress=unmatchedSuppression
  --suppress=passedByValue
  --suppress=syntaxError
  --inconclusive
)

target_compile_features(project_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})
add_library(autoDMP::project_options INTERFACE IMPORTED)
add_library(autoDMP::project_warnings INTERFACE IMPORTED)

# configure files based on CMake configuration options
add_subdirectory(configured_files)
# target_include_directories(project_options INTERFACE "${CMAKE_BINARY_DIR}")

# Adding src directories:
add_subdirectory(src)
add_subdirectory(apps)

# Adding the tests:
# option(ENABLE_TESTING "Enable the tests" ${PROJECT_IS_TOP_LEVEL})
# if(ENABLE_TESTING)
#   enable_testing()
#   message(AUTHOR_WARNING "Building Tests. Be sure to check out test/constexpr_tests.cpp for constexpr testing")
#   add_subdirectory(test)
# endif()

# option(ENABLE_FUZZING "Enable the fuzz tests" OFF)
# if(ENABLE_FUZZING)
#   message(AUTHOR_WARNING "Building Fuzz Tests, using fuzzing sanitizer https://www.llvm.org/docs/LibFuzzer.html")
#   add_subdirectory(fuzz_test)
# endif()

if(CMAKE_SKIP_INSTALL_RULES)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  return()
endif()

# Add other targets that you want installed here, be default we just package the one executable
# we know we want to ship
package_project(TARGETS intro project_options project_warnings)

# Experience shows that explicit package naming can help make it easier to sort
# out potential ABI related issues before they start, while helping you
# track a build to a specific GIT SHA
set(CPACK_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}-${GIT_SHORT_SHA}-${CMAKE_SYSTEM_NAME}-${CMAKE_BUILD_TYPE}-${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}"
)

include(CPack)
