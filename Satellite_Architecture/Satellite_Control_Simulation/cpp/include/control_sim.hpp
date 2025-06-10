#pragma once
#include <array>
#include <vector>
#include "euler_dynamics.hpp"

using Vec3 = std::array<double,3>;
using Vec4 = std::array<double,4>;

/// Attitude state: body rates + quaternion
struct State {
    Vec3 omega;
    Vec4 quat;
};

/// History containers
struct History {
    std::vector<State> states;
    std::vector<Vec3>  control_torques;
};

/**
 * Fixed‐step RK4 simulation of attitude control.
 *
 * @param t0     start time [s]
 * @param tf     end time   [s]
 * @param dt     step size  [s]
 * @param x0     initial state
 * @param k_bdot true → bang‐bang B‐dot controller
 * @returns history of states & torques
 */
History simulateControl(double t0, double tf, double dt,
                        const State &x0, bool k_bdot);
