// include/adcs_hal.h
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

    /// Which control law should run
    typedef enum {
        CTRL_MODE_B_DOT = 0,  ///< classic B·dot bang–bang
        CTRL_MODE_WX_B   = 1   ///< ω×B law (despin)
    } ControlMode;

    /// Sensor readings passed into your control algorithm each step
    typedef struct {
        double wx, wy, wz;   ///< body‐rates (rad/s)
        double Bx, By, Bz;   ///< magnetic field in body frame (T)
    } SensorReadings;

    /// Dipole commands your algorithm outputs each step
    typedef struct {
        double mx, my, mz;   ///< magnetic dipole moment (A·m²)
    } ActuatorCommands;

    /// Initialize the HAL (called once at sim startup)
    /// @param sample_interval_s control‐step interval, in seconds
    void HAL_init(double sample_interval_s);

    /// The core control entry point: your embedded code implements this.
    /// @param sens  pointer to current sensor readings
    /// @param act   pointer to output actuator commands
    void HAL_compute_control(const SensorReadings* sens,
                             ActuatorCommands*   act);

    /// Block/delay for @p seconds (simulate compute time or RTOS delay)
    void HAL_delay(double seconds);

    /// Switch between your two control modes at runtime
    /// @param mode  one of CTRL_MODE_B_DOT or CTRL_MODE_WX_B
    void HAL_set_control_mode(ControlMode mode);

#ifdef __cplusplus
}
#endif
