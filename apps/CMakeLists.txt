add_executable(autoDMP "./src/main.cpp")
target_link_libraries(autoDMP PRIVATE project_options project_warnings)

# includes (private includes require absolute paths)
set(autoDMP_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  CACHE STRING "")
target_include_directories(autoDMP PRIVATE ${autoDMP_INCLUDE_DIR})

# # find dependencies (disabled for now since no dependencies need to be found)
# set(autoDMP_DEPENDENCIES_CONFIGURED
#     fmt
#     CACHE STRING "")
# foreach (dep ${autoDMP_DEPENDENCIES_CONFIGURED})
#   find_package(${dep} CONFIG REQUIRED)
# endforeach()

# link dependencies
set(autoDMP_LINKED_LIBRARIES util) # cam gui imcap improc)
target_link_system_libraries(autoDMP PRIVATE ${autoDMP_LINKED_LIBRARIES})

# tests
if(FEATURE_TESTS)
  add_subdirectory("./test")
endif()
