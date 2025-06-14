# src/CMakeLists.txt

# Gather all simulation source files
file(GLOB SIM_SRCS
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

# Define the simulator executable
add_executable(SatelliteSim
        ${SIM_SRCS}
)

# Link against matio (for .mat IGRF coefficient reading)
target_link_libraries(SatelliteSim
        PRIVATE
        PkgConfig::matio
)

# Expose the chosen output directory to the C++ code:
# This defines a macro SIM_OUTPUT_DIR which expands to your SAT_SIM_OUTPUT_DIR path.
target_compile_definitions(SatelliteSim
        PRIVATE
        SIM_OUTPUT_DIR="${SAT_SIM_OUTPUT_DIR}"
)

# (Optional) If you need filesystem, etc.
target_link_libraries(SatelliteSim
        PRIVATE
        stdc++fs    # on older toolchains; remove if not needed
)
