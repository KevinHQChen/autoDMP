# include guard
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
  return()
endif()

set(ENABLE_CLANG_TIDY OFF)
set(ENABLE_CPPCHECK OFF)
set(ENABLE_SANITIZER_ADDRESS OFF)
set(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR OFF)
set(ENABLE_COVERAGE OFF)

# building the tests (check test/constexpr_test.cpp for constexpr testing)
# enable sanitizers and clang-tidy if running the tests
option(FEATURE_TESTS "Enable the tests" OFF)
if(FEATURE_TESTS)
  list(APPEND VCPKG_MANIFEST_FEATURES "tests") # activate the tests feature in the manifest (vcpkg.json)
endif()
if(FEATURE_TESTS)
  set(ENABLE_CLANG_TIDY "ENABLE_CLANG_TIDY")
  set(ENABLE_CPPCHECK "ENABLE_CPPCHECK")
  set(ENABLE_COVERAGE "ENABLE_COVERAGE")
  set(ENABLE_SANITIZER_ADDRESS "ENABLE_SANITIZER_ADDRESS")
  set(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR "ENABLE_SANITIZER_UNDEFINED_BEHAVIOR")
endif()

# building the docs
option(FEATURE_DOCS "Enable the docs" OFF)
if(FEATURE_DOCS)
  set(ENABLE_DOXYGEN "ENABLE_DOXYGEN")
else()
  set(ENABLE_DOXYGEN OFF)
endif()

# fuzz tests (using fuzzing sanitizer libFuzzer)
option(FEATURE_FUZZ_TESTS "Enable the fuzz tests" OFF)
if(FEATURE_FUZZ_TESTS)
  add_subdirectory(fuzz_test)
endif()

# Initialize project_options variable related to this project
# This overwrites `project_options` and sets `project_warnings`
# uncomment to enable the options. Some of them accept one or more inputs:
# project_options(
#   ENABLE_CACHE
#   ${ENABLE_CPPCHECK}
#   ${ENABLE_CLANG_TIDY}
#   # ENABLE_INTERPROCEDURAL_OPTIMIZATION
#   # ENABLE_NATIVE_OPTIMIZATION
#   ${ENABLE_DOXYGEN}
#   ${ENABLE_COVERAGE}
#   ${ENABLE_SANITIZER_ADDRESS}
#   # ENABLE_SANITIZER_UNDEFINED_BEHAVIOR
#   # ENABLE_SANITIZER_LEAK
#   # ENABLE_SANITIZER_THREAD
#   # ENABLE_SANITIZER_MEMORY
#   # ENABLE_PCH
#   # PCH_HEADERS
#   # WARNINGS_AS_ERRORS
#   # ENABLE_INCLUDE_WHAT_YOU_USE
#   # ENABLE_USER_LINKER
#   # ENABLE_BUILD_WITH_TIME_TRACE
#   # ENABLE_UNITY
# )
