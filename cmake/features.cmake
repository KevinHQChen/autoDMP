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

# Enable CCache if available
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

# building the tests (check test/constexpr_test.cpp for constexpr testing)
# enable sanitizers and clang-tidy if running the tests
option(FEATURE_TESTS "Enable the tests" OFF)
if(FEATURE_TESTS)
  list(APPEND VCPKG_MANIFEST_FEATURES "tests") # activate the tests feature in the manifest (vcpkg.json)
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
