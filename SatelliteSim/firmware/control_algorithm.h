#pragma once

#include "adcs_hal.h"  // brings in ControlMode, SensorReadings, ActuatorCommands

#ifdef __cplusplus
extern "C" {
#endif

    /// Configure which control law runs (B-dot or ω×B)
    void control_algorithm_set_mode(ControlMode mode);

    /// Main control entry point: reads sensors, writes actuator commands
    void control_algorithm(const SensorReadings* sens,
                           ActuatorCommands*   act);

#ifdef __cplusplus
}
#endif
