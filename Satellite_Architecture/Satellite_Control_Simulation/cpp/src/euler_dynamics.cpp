#include "euler_dynamics.hpp"
#include <cmath>

Vec7 eulerDynamics(double /*t*/, const Vec7 &y, const Vec3 &L) {
    // Inertia constants (same as MATLAB):
    constexpr double mass = 10.0;
    constexpr double h = 0.34, w = 0.2, d = 0.1;
    constexpr double Jx = mass/12*(w*w + d*d);
    constexpr double Jy = mass/12*(h*h + d*d);
    constexpr double Jz = mass/12*(h*h + w*w);

    double w1 = y[0], w2 = y[1], w3 = y[2];
    double q1 = y[3], q2 = y[4], q3 = y[5], q4 = y[6];

    // Angular acceleration: ω̇ = I⁻¹ (L - ω×(Iω))
    double wd1 = (L[0] - (Jy - Jz)*w2*w3) / Jx;
    double wd2 = (L[1] - (Jz - Jx)*w3*w1) / Jy;
    double wd3 = (L[2] - (Jx - Jy)*w1*w2) / Jz;

    // Quaternion kinematics: q̇ = 0.5 * Ω(ω) * q
    Vec7 dydt;
    dydt[0] = wd1;
    dydt[1] = wd2;
    dydt[2] = wd3;
    dydt[3] =  0.5*(  q4*w1 - q3*w2 + q2*w3);
    dydt[4] =  0.5*(  q3*w1 + q4*w2 - q1*w3);
    dydt[5] =  0.5*( -q2*w1 + q1*w2 + q4*w3);
    dydt[6] =  0.5*( -q1*w1 - q2*w2 - q3*w3);

    return dydt;
}