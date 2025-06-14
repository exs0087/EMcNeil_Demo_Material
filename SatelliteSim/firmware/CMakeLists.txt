# firmware/CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(control_firmware C)

# Expose include directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Build your embedded control logic as a static library
add_library(control_algorithm STATIC
        control_algorithm.c
)

# Install headers if desired
install(FILES control_algorithm.h DESTINATION include)
