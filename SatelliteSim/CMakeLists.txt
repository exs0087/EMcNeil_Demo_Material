cmake_minimum_required(VERSION 3.15)
project(SatelliteSim LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------
# Include directories
include_directories(
        ${CMAKE_SOURCE_DIR}/include   # your headers (e.g. iss_orbit.h)
        ${CMAKE_SOURCE_DIR}/firmware
)

# ------------------------------------------------------------
# control_algorithm static lib
add_library(control_algorithm
        firmware/control_algorithm.cpp
)

# ------------------------------------------------------------
# Core simulation library
add_library(sim_core
        src/adcs_hal_sim.cpp
        src/eulerseqns2.cpp
        src/igrf.cpp
        src/igrf_loader.cpp
        src/iss_orbit.cpp        # only the .cpp, not the .h!
        src/magnetic_field.cpp
        src/ned2eci.cpp
        src/ode4.cpp
        src/vendor_timer_stub.cpp
)
target_link_libraries(sim_core PUBLIC control_algorithm)

# ------------------------------------------------------------
# Simulator executable
add_executable(simulator
        src/main.cpp
)
target_link_libraries(simulator PRIVATE sim_core)

# ------------------------------------------------------------
# Test executable
add_executable(test_bfield
        src/test_bfield.cpp
)
target_link_libraries(test_bfield PRIVATE sim_core)

# ------------------------------------------------------------
# Allow overriding T_FINAL at configure time
if (DEFINED T_FINAL)
    target_compile_definitions(simulator PRIVATE T_FINAL=${T_FINAL})
endif()

# ------------------------------------------------------------
# Install rules (optional)
install(TARGETS simulator DESTINATION bin)
