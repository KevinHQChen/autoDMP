# good cmake resource: https://cliutils.gitlab.io/modern-cmake/, https://alexreinking.com/blog/how-to-use-cmake-without-the-agonizing-pain-part-1.html
# include guard
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
  return()
endif()

# Enable CCache if available (https://crascit.com/2016/04/09/using-ccache-with-cmake/)
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

# enable building the tests (check test/constexpr_test.cpp for constexpr testing)
option(BUILD_TESTS "Enable building tests" OFF)
option(ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
option(ENABLE_CPPCHECK "Enable cppcheck" OFF)
option(ENABLE_COVERAGE "Enable coverage" OFF)
option(ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" OFF)
option(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR "Enable undefined behavior sanitizer" OFF)
if(BUILD_TESTS)
  list(APPEND VCPKG_MANIFEST_FEATURES "tests") # activate the tests feature in the manifest (vcpkg.json)
  set(ENABLE_CLANG_TIDY "ENABLE_CLANG_TIDY")
  set(ENABLE_CPPCHECK "ENABLE_CPPCHECK")
  set(ENABLE_COVERAGE "ENABLE_COVERAGE")
  set(ENABLE_SANITIZER_ADDRESS "ENABLE_SANITIZER_ADDRESS")
  set(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR "ENABLE_SANITIZER_UNDEFINED_BEHAVIOR")
  enable_testing()
endif()

# enable linting and static analysis if building tests
if(ENABLE_CLANG_TIDY)
  find_program(CLANG_TIDY_PROGRAM clang-tidy)
  if(CLANG_TIDY_PROGRAM)
    set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_PROGRAM}")
  endif()
endif()
if(ENABLE_CPPCHECK)
  find_program(CPPCHECK_PROGRAM cppcheck)
  if(CPPCHECK_PROGRAM)
    set(CMAKE_CXX_CPPCHECK "${CPPCHECK_PROGRAM}" "--enable=all")
  endif()
endif()

# TODO enable coverage if building tests

# TODO enable ASAN and UBSAN if building tests (https://stackoverflow.com/questions/44320465/whats-the-proper-way-to-enable-addresssanitizer-in-cmake-that-works-in-xcode)

# enable building the docs
# (first run `doxygen -g Doxyfile` to generate a default doxygen config file,
# set RECURSIVE to YES and add the relevant input files to the INPUT field)
option(BUILD_DOCS "Enable the docs" OFF)
if(BUILD_DOCS)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/docs/Doxyfile)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    add_custom_target(
      docs
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen"
      VERBATIM)
  else()
    message(WARNING "Doxygen not found, unable to generate documentation.")
  endif()
endif()

# TODO fuzz tests (using fuzzing sanitizer libFuzzer)
option(FEATURE_FUZZ_TESTS "Enable the fuzz tests" OFF)
if(FEATURE_FUZZ_TESTS)
  add_subdirectory(fuzz_test)
endif()

# cool features to add for fun:
# "WARNINGS_AS_ERRORS\;OFF\;ON\;Treat warnings as Errors"
# "ENABLE_SANITIZER_ADDRESS\;OFF\;${SUPPORTS_ASAN}\;Make memory errors into hard runtime errors (windows/linux/macos)"
# "ENABLE_SANITIZER_UNDEFINED_BEHAVIOR\;OFF\;${SUPPORTS_UBSAN}\;Make certain types (numeric mostly) of undefined behavior into runtime errors"
# "ENABLE_INTERPROCEDURAL_OPTIMIZATION\;OFF\;OFF\;Enable whole-program optimization (e.g. LTO)"
# "ENABLE_NATIVE_OPTIMIZATION\;OFF\;OFF\;Enable the optimizations specific to the build machine (e.g. SSE4_1, AVX2, etc.)."
# "DISABLE_EXCEPTIONS\;OFF\;OFF\;Disable Exceptions (no-exceptions and no-unwind-tables flag)"
# "DISABLE_RTTI\;OFF\;OFF\;Disable RTTI (no-rtti flag)"
# "ENABLE_INCLUDE_WHAT_YOU_USE\;OFF\;OFF\;Enable include-what-you-use analysis during compilation"
# "ENABLE_PCH\;OFF\;OFF\;Enable pre-compiled-headers support"
# "ENABLE_BUILD_WITH_TIME_TRACE\;OFF\;OFF\;Generates report of where compile-time is spent"
# "ENABLE_UNITY\;OFF\;OFF\;Merge C++ files into larger C++ files, can speed up compilation sometimes"
# "ENABLE_SANITIZER_LEAK\;OFF\;OFF\;Make memory leaks into hard runtime errors"
# "ENABLE_SANITIZER_THREAD\;OFF\;OFF\;Make thread race conditions into hard runtime errors"
# "ENABLE_SANITIZER_MEMORY\;OFF\;OFF\;Make other memory errors into runtime errors")
