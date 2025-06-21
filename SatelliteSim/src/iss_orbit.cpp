// src/iss_orbit.cpp
// Updated to include geodetic latitude, longitude, and altitude computation

#include "iss_orbit.h"
#include "Vec3.h"
#include <cmath>

// Earth and orbit constants
static constexpr double mu          = 398600.0;            // km^3/s^2
static constexpr double Re_km       = 6371.0;               // km
static constexpr double alt_km_nom  = 400.0;                // km (ISS altitude)
static constexpr double inclination = 52.0 * M_PI/180.0;    // rad
static constexpr double theta0      = -M_PI/2;              // rad, initial true anomaly

// Derived orbital parameters
static const double r_orbit   = Re_km + alt_km_nom;                       // km
static const double T_orbit   = 2*M_PI * std::pow(r_orbit,1.5) / std::sqrt(mu); // s
static const double theta_dot = 2*M_PI / T_orbit;                         // rad/s

//------------------------------------------------------------------------------
// Compute ECEF position (km) of ISS at time t, then convert to geodetic
// latitude (rad), longitude (rad), and altitude above ellipsoid (km).
//------------------------------------------------------------------------------
void issOrbitGeod(double t, double& lat, double& lon, double& alt_km) {
    // 1) Inertial position in orbital plane
    double theta = theta0 + theta_dot * t;
    double x_orb = r_orbit * std::cos(theta);
    double y_orb = r_orbit * std::sin(theta);

    // 2) Rotate by inclination about X-axis into ECI
    double x_eci = x_orb;
    double y_eci = y_orb * std::cos(inclination);
    double z_eci = y_orb * std::sin(inclination);

    // 3) Rotate ECI -> ECEF (Earth rotation)
    double gst    = we * t;
    double ce     = std::cos(gst);
    double se     = std::sin(gst);
    double x_ecef =  ce * x_eci + se * y_eci;
    double y_ecef = -se * x_eci + ce * y_eci;
    double z_ecef =  z_eci;

    // 4) Compute geocentric radius
    double r = std::hypot(x_ecef, y_ecef, z_ecef);

    // 5) Longitude
    lon = std::atan2(y_ecef, x_ecef);

    // 6) Geodetic latitude & altitude via WGS-84 ellipsoid
    //    (a = equatorial radius, f = flattening)
    constexpr double a = 6378.137;            // km
    constexpr double f = 1.0/298.257223563;
    double b = a * (1 - f);

    // Compute initial approximation
    double p = std::hypot(x_ecef, y_ecef);
    double th = std::atan2(a * z_ecef, b * p);
    double sin_th = std::sin(th), cos_th = std::cos(th);
    lat = std::atan2(z_ecef + (f*(2 - f))*b*std::pow(sin_th,3),
                     p - (f*(2 - f))*a*std::pow(cos_th,3));

    // Radius of curvature in prime vertical
    double N = a / std::sqrt(1 - (2*f - f*f)*std::sin(lat)*std::sin(lat));
    alt_km = p / std::cos(lat) - N;
}

// Alias to match existing API (returns only ECI vector)
Vec3 issOrbitECI(double t) {
    // We still want the raw ECI coords for other uses:
    double theta = theta0 + theta_dot * t;
    double x_orb = r_orbit * std::cos(theta);
    double y_orb = r_orbit * std::sin(theta);

    double x = x_orb;
    double y = y_orb * std::cos(inclination);
    double z = y_orb * std::sin(inclination);

    double gst = we * t;
    double c = std::cos(gst), s = std::sin(gst);
    return Vec3{
        c*x + s*y,
       -s*x + c*y,
        z
    };
}
