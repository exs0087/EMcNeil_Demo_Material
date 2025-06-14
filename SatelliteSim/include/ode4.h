#pragma once

#include <vector>
#include <functional>
#include "Vec7.h"

/// Simple fixed-step RK4 integrator.
/// @param f      Function f(t, y) returning dy/dt.
/// @param t      Vector of time points (must be strictly increasing).
/// @param y0     Initial state at t[0].
/// @returns      Vector Y of size t.size(), with Y[i] ≈ y(t[i]).
std::vector<Vec7> ode4(
    const std::function<Vec7(double, const Vec7&)>& f,
    const std::vector<double>& t,
    const Vec7& y0
);
