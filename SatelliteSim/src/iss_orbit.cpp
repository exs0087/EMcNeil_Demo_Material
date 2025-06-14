// src/iss_orbit.cpp

#include "iss_orbit.h"
#include "Vec3.h"
#include <cmath>

// Port of your MATLAB issorb.m, but all in meters & SI units:
Vec3 issOrbitECI(double t) {
    constexpr double pi       = 3.141592653589793;
    constexpr double incl     = 52.0 * pi / 180.0;   // inclination [rad]
    constexpr double theta0   = -pi/2.0;             // initial anomaly [rad]
    constexpr double alt_m    = 400e3;               // altitude [m]
    constexpr double Re_m     = 6371e3;              // Earth radius [m]
    constexpr double mu       = 3.986004418e14;      // GM [m^3/s^2]

    double r     = Re_m + alt_m;
    double theta = theta0 + std::sqrt(mu/(r*r*r)) * t;  // rad/s * s

    double c0 = std::cos(incl), s0 = std::sin(incl);
    double ct = std::cos(theta), st = std::sin(theta);

    return Vec3{
        r * ct *  c0,
        r * st *  c0,
        r * st *  s0
    };
}
