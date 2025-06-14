#pragma once

#include "Vec7.h"
#include "adcs_hal.h"    // for SensorReadings, ActuatorCommands
#include <vector>

/**
 * Computes the time‐derivative for the 7‐element state vector:
 *   [ ω₁, ω₂, ω₃, q₁, q₂, q₃, q₄ ]
 * Implements rigid‐body Euler + quaternion kinematics and calls into
 * the on‐board controller via the HAL.
 */
Vec7 eulerseqns2(double t, const Vec7& y);

/// Get the history of actuator commands (one entry per top‐level time‐step).
const std::vector<ActuatorCommands>& getActHistory();

/// Clear out any previously recorded actuator history (call before each run).
void resetActHistory();
