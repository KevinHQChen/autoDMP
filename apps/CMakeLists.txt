add_executable(autoDMP "./src/main.cpp")
# target_link_libraries(autoDMP PRIVATE project_options project_warnings)

# includes (private includes require absolute paths)
set(autoDMP_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  CACHE STRING "")
target_include_directories(autoDMP PRIVATE ${autoDMP_INCLUDE_DIR})

# link dependencies
target_link_options(autoDMP PUBLIC "LINKER:-export-dynamic")
target_link_system_libraries(autoDMP PRIVATE
  util
  cam
  pump
  gui
  imcap
  improc
  ctrl
  imgui_bundle
)

# tests
if(FEATURE_TESTS)
  add_subdirectory("./test")
endif()
