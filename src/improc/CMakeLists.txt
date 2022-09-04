# Image processing library for autoDMP
set(HEADER_LIST "../../include/improc/improc.hpp")

# Make an automatic library - will be static or dynamic based on user setting
add_library(improc improc.cpp ${HEADER_LIST})     # target is created


target_link_libraries(improc
    PUBLIC util)    # link each library where the target is created

# We need this directory, and users of our library will need it too
target_include_directories(improc
    PUBLIC ../../include)

# Image capture library for autoDMP
set(HEADER_LIST "../../include/improc/imcap.hpp")

add_library(imcap imcap.cpp ${HEADER_LIST})     # target is created

target_link_libraries(imcap
    PUBLIC cam
    PUBLIC util)    # link each library where the target is created

# We need this directory, and users of our library will need it too
target_include_directories(imcap
    PUBLIC ../../include)
