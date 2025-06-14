// firmware/control_algorithm.h
#pragma once
#include "../include/adcs_hal.h"  // for ControlMode, SensorReadings, ActuatorCommands

#ifdef __cplusplus
extern "C" {
#endif

    /// Set the active control mode.  (Called from sim stub or CLI.)
    void control_algorithm_set_mode(ControlMode mode);

    /// Your on‐board attitude control algorithm entry point.
    void control_algorithm(const SensorReadings* sens,
                           ActuatorCommands*   act);

#ifdef __cplusplus
}
#endif
