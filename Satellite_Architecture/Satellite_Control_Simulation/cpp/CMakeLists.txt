cmake_minimum_required(VERSION 3.15)
project(satellite_control LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Fetch pybind11 (for your Python bindings)
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.11.1
)
FetchContent_MakeAvailable(pybind11)

# Tell CMake where to find your own headers:
include_directories(${PROJECT_SOURCE_DIR}/include)

# -----------------------------------
# 1) Build the Python extension module
# -----------------------------------
pybind11_add_module(satellite_kernels
    src/euler_dynamics.cpp
    src/control_sim.cpp
    src/magnetic_field.cpp
    src/bindings.cpp
)

# -----------------------------------
# 2) Build the standalone C++ simulator
# -----------------------------------
# List all .cpp files needed for the despin simulation:
set(DESPIN_SRCS
    src/main.cpp
    src/ode4m.cpp
    src/eulerseqns2.cpp
    src/iss_orbit.cpp
    src/igrf.cpp
    src/IGRFLoader.cpp
    src/legendre.cpp
    src/ned2eci.cpp
)

add_executable(despin_sim ${DESPIN_SRCS})

# If any external libraries are needed, link them here:
# target_link_libraries(despin_sim PRIVATE <other_libraries>)

# -----------------------------------
# 3) Install or copy datfiles folder
# -----------------------------------
# The IGRFLoader expects a folder named "datfiles/" next to the executable.
# We’ll copy it into the build output directory automatically.
add_custom_command(TARGET despin_sim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${PROJECT_SOURCE_DIR}/datfiles
    $<TARGET_FILE_DIR:despin_sim>/datfiles
)

# -----------------------------------
# 4) (Optional) Install step
# -----------------------------------
install(TARGETS despin_sim RUNTIME DESTINATION bin)
install(TARGETS satellite_kernels LIBRARY DESTINATION python)

