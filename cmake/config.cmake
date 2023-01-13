# TODO Don't use dynamic project options, since we don't even know why the default options are the way they are, we should play around with them manually to see which ones are actually useful for us

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
