add_library(imcap "./src/imcap/imcap.cpp")
# target_link_libraries(imcap PRIVATE project_options project_warnings)

# includes
set(imcap_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  CACHE STRING "")
target_include_directories(imcap PUBLIC
    "$<BUILD_INTERFACE:${imcap_INCLUDE_DIR}>"
    "$<INSTALL_INTERFACE:./${CMAKE_INSTALL_INCLUDEDIR}>")

target_link_system_libraries(imcap PRIVATE
  util
  cam)

# tests
if(FEATURE_TESTS)
  add_subdirectory("./test")
endif()
