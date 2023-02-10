add_library(util "./src/util/util.cpp")
# target_link_libraries(util PRIVATE project_options project_warnings)

# includes
set(util_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  CACHE STRING "")
target_include_directories(util PUBLIC
  "$<BUILD_INTERFACE:${util_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:./${CMAKE_INSTALL_INCLUDEDIR}>")

# find dependencies
set(util_DEPENDENCIES_CONFIGURED
    fmt
    spdlog
    CLI11
    toml11
    tsl-ordered-map
    OpenCV
    CACHE STRING "")
foreach (dep ${util_DEPENDENCIES_CONFIGURED})
  find_package(${dep} REQUIRED) # certain packages like boost may need different find_package calls, check the output of run_vcpkg
  if (DEFINED ${dep}_INCLUDE_DIRS)
    target_include_directories(util PUBLIC ${${dep}_INCLUDE_DIRS})
  endif ()
endforeach ()

# link dependencies
target_link_system_libraries(util PRIVATE
    # unofficial::nativefiledialog::nfd
    fmt::fmt fmt::fmt-header-only
    spdlog::spdlog spdlog::spdlog_header_only
    CLI11::CLI11
    toml11::toml11
    tsl::ordered_map
    ${OpenCV_LIBS}
    # Boost::program_options
)

# tests
if(FEATURE_TESTS)
  add_subdirectory("./test")
endif()
