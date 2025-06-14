#include "adcs_hal.h"
#include "control_algorithm.h"   // for control_algorithm_set_mode()
#include <chrono>
#include <thread>

static double g_sample_interval_s = 0.0;

extern "C" {

    void HAL_init(double sample_interval_s) {
        g_sample_interval_s = sample_interval_s;
        // Could reset any sim‐specific buffers here
    }

    void HAL_compute_control(const SensorReadings* sens,
                             ActuatorCommands*   act) {
        // Call into your real embedded algorithm
        control_algorithm(sens, act);
    }

    void HAL_delay(double seconds) {
        // Simulate compute/RTOS delay
        std::this_thread::sleep_for(
            std::chrono::duration<double>(seconds)
        );
    }

    void HAL_set_control_mode(ControlMode mode) {
        // Route the sim’s mode switch into your firmware stub
        control_algorithm_set_mode(mode);
    }

}  // extern "C"
