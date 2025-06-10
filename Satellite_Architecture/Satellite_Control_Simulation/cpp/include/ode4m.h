// ode4m.h
#ifndef ODE4M_H
#define ODE4M_H

#include <vector>
#include <functional>
#include <array>

// ode4m: 4th-order fixed-step integrator with magnetorquer controls
// Template: state vector size must be 7 (3x rates, 4x quaternion)
// odefun: user-provided dynamics function: f(t, y) -> dy/dt (size 7)
// tspan: vector of times [t0, t1, ..., tN]
// y0: initial state (length 7)
// cmode: control mode (1 = B-dot, 2 = omega x B)
// Returns: solution Y (size Nx7), each row corresponds to time in tspan

std::vector<std::array<double,7>> ode4m(
    const std::function<std::array<double,7>(double, const std::array<double,7>&)> &odefun,
    const std::vector<double> &tspan,
    const std::array<double,7> &y0,
    int cmode
);

#endif // ODE4M_H