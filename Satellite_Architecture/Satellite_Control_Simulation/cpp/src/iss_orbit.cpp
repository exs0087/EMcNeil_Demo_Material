#include "iss_orbit.h"
#include <cmath>

void issOrbit(double t, double &latitude, double &longitude, double &altitude) {
    // Convert degrees → radians
    constexpr double deg2rad = M_PI / 180.0;
    // Inclination = 52° → radians
    const double alpha = 52.0 * deg2rad;
    // Initial orbital argument (theta at t=0)
    double theta = -M_PI / 2.0;
    // Fixed altitude above Earth [km]
    const double alt_fixed = 400.0;

    // Standard gravitational parameter and Earth radius
    const double mu = 398600.0;  // km^3/s^2
    const double Re = 6371.0;    // km

    // Orbital radius
    const double r = Re + alt_fixed; // km

    // Circular orbital period: T = 2*pi*sqrt(r^3/μ)
    const double T_orbit = 2.0*M_PI * std::sqrt(r*r*r / mu);
    // Angular rate [rad/s]
    const double theta_dot = 2.0*M_PI / T_orbit;

    // Advance theta by t
    theta += theta_dot * t;

    // ECI coordinates of a circular, inclined orbit:
    // x = r * cos(theta) * cos(alpha)
    // y = r * sin(theta) * cos(alpha)
    // z = r * sin(theta) * sin(alpha)
    double x = r * std::cos(theta) * std::cos(alpha);
    double y = r * std::sin(theta) * std::cos(alpha);
    double z = r * std::sin(theta) * std::sin(alpha);

    // Project into equatorial plane
    double p = std::sqrt(x*x + y*y);

    // Geocentric latitude = atan2(z, p)
    latitude = std::atan2(z, p);

    // Longitude in [0, 2*pi)
    double lon = std::atan2(y, x);
    if (lon < 0.0) {
        lon += 2.0 * M_PI;
    }
    longitude = lon;

    // Altitude = constant 400 km
    altitude = alt_fixed;
}
