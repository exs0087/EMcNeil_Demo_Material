// src/adcs_hal_sim.cpp

#include "adcs_hal.h"
#include "control_algorithm.h"   // for control_algorithm_set_mode()
#include <chrono>
#include <thread>

// Simulation-wide sample interval (in seconds)
static double g_sample_interval_s = 0.0;

extern "C" {

    // Initialize the hardware abstraction layer with a given sample interval [seconds].
    void HAL_init(double sample_interval_s) {
        g_sample_interval_s = sample_interval_s;
        // If you want to reset any sim-specific state, do so here.
    }

    // Compute actuator commands given current sensor readings (all SI units!).
    void HAL_compute_control(const SensorReadings* sens,
                             ActuatorCommands*   act) {
        // Forward directly to the embedded/flight algorithm.
        // Sensor readings: ω [rad/s], B [Tesla], time [s]
        control_algorithm(sens, act);
    }

    // Simulate a hardware delay (e.g., control loop period) [seconds].
    void HAL_delay(double seconds) {
        std::this_thread::sleep_for(
            std::chrono::duration<double>(seconds)
        );
    }

    // Set the control mode (e.g., B-dot vs ω×B)
    void HAL_set_control_mode(ControlMode mode) {
        // Route the sim's mode switch into the firmware stub
        control_algorithm_set_mode(mode);
    }

}  // extern "C"
