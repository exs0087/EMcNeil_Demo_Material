cmake_minimum_required(VERSION 3.10)
project(SatelliteSim LANGUAGES C CXX)

# 1) Require C++17
set(CMAKE_CXX_STANDARD      17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#
# 2) Build your on-board control algorithm as a static lib
#
add_library(control_algorithm STATIC
        firmware/control_algorithm.cpp
)
target_include_directories(control_algorithm PUBLIC
        ${PROJECT_SOURCE_DIR}/include    # for adcs_hal.h, control_algorithm.h
        ${PROJECT_SOURCE_DIR}/firmware   # for control_algorithm.h
)

#
# 3) Build the simulator executable
#
add_executable(simulator
        src/adcs_hal_sim.cpp
        src/eulerseqns2.cpp
        src/igrf.cpp
        src/igrf_loader.cpp
        src/iss_orbit.cpp
        src/magnetic_field.cpp
        src/main.cpp
        src/ned2eci.cpp
        src/ode4.cpp
)
target_include_directories(simulator PRIVATE
        ${PROJECT_SOURCE_DIR}/include    # for all our public headers
)
target_link_libraries(simulator PRIVATE
        control_algorithm                # hook in your firmware lib
)

#
# 4) Locate and link libmatio (no pkg-config)
#
find_library(MATIO_LIB matio
        HINTS /usr/lib /usr/local/lib
)
find_path(MATIO_INC matio.h
        HINTS /usr/include /usr/local/include
)
if(NOT MATIO_LIB OR NOT MATIO_INC)
    message(FATAL_ERROR "libmatio not found; install libmatio-dev and retry")
endif()

target_include_directories(simulator PRIVATE ${MATIO_INC})
target_link_libraries(simulator      PRIVATE ${MATIO_LIB})

# On Unix we also need the math library
target_link_libraries(simulator      PRIVATE m)
