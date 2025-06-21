// src/ode4.cpp

#include "ode4.h"
#include <cmath>

// Classic RK4 integrator with built‐in quaternion renormalization
std::vector<Vec7> ode4(
    const std::function<Vec7(double, const Vec7&)>& f,
    const std::vector<double>& t,
    const Vec7& y0
) {
    const size_t N = t.size();
    std::vector<Vec7> Y(N);
    Y[0] = y0;

    for (size_t i = 1; i < N; ++i) {
        double h      = t[i]     - t[i - 1];
        double t_prev = t[i - 1];
        const Vec7& y_prev = Y[i - 1];

        // Four slopes
        Vec7 k1 = f(t_prev,              y_prev);
        Vec7 k2 = f(t_prev + 0.5 * h,     y_prev + k1 * (0.5 * h));
        Vec7 k3 = f(t_prev + 0.5 * h,     y_prev + k2 * (0.5 * h));
        Vec7 k4 = f(t_prev +       h,     y_prev + k3 * h);

        // Advance
        Vec7 y = y_prev + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);

        // Renormalize quaternion part (indices 3–6)
        double q1 = y[3];
        double q2 = y[4];
        double q3 = y[5];
        double q4 = y[6];
        double qnorm = std::sqrt(
            q1*q1 +
            q2*q2 +
            q3*q3 +
            q4*q4
        );

        if (qnorm > 0.0) {
            y[3] = q1 / qnorm;
            y[4] = q2 / qnorm;
            y[5] = q3 / qnorm;
            y[6] = q4 / qnorm;
        } else {
            // Should never happen unless you hit NaN earlier;
            // reset to identity quaternion
            y[3] = 0.0;
            y[4] = 0.0;
            y[5] = 0.0;
            y[6] = 1.0;
        }

        Y[i] = y;
    }

    return Y;
}
