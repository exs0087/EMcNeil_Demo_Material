// src/magnetic_field.cpp

#include "magnetic_field.h"
#include "igrf.h"
#include "ned2eci.h"
#include "Vec3.h"
#include "Vec4.h"
#include "iss_orbit.h"
#include <cmath>

// Earth’s sidereal rotation rate [rad/s]
static constexpr double omegaEarth = 7.2921150e-5;

// Base date: June 12 2025 as decimal year
static constexpr double BASE_YEAR = 2025.0 + (163.0 - 1.0)/365.25;

// Rotate a vector from ECI ➔ ECEF (positions in km)
static Vec3 eciToEcef(const Vec3& v, double t_sec) {
    double θ = omegaEarth * t_sec;
    double c = std::cos(θ), s = std::sin(θ);
    return Vec3{
        + c*v.x + s*v.y,
        - s*v.x + c*v.y,
          v.z
    };
}

// Rotate a vector from ECEF ➔ ECI
static Vec3 ecefToEci(const Vec3& v, double t_sec) {
    double θ = omegaEarth * t_sec;
    double c = std::cos(θ), s = std::sin(θ);
    return Vec3{
          c*v.x - s*v.y,
          s*v.x + c*v.y,
          v.z
    };
}

// Build body‐to‐ECI DCM from quaternion (w last), then transpose:
static Vec3 rotateEciToBody(const Vec3& v, const Vec4& q) {
    double w=q.w, x=q.x, y=q.y, z=q.z;
    double C11 =  1 - 2*(y*y + z*z);
    double C12 =      2*(x*y + z*w);
    double C13 =      2*(x*z - y*w);
    double C21 =      2*(x*y - z*w);
    double C22 =  1 - 2*(x*x + z*z);
    double C23 =      2*(y*z + x*w);
    double C31 =      2*(x*z + y*w);
    double C32 =      2*(y*z - x*w);
    double C33 =  1 - 2*(x*x + y*y);

    // transpose to go ECI ➔ body
    return Vec3{
        C11*v.x + C21*v.y + C31*v.z,
        C12*v.x + C22*v.y + C32*v.z,
        C13*v.x + C23*v.y + C33*v.z
    };
}

Vec3 getBodyFieldAt(double t_sec, const Vec4& q) {
    // 1) Compute decimal year
    double year = BASE_YEAR + t_sec/(86400.0*365.25);

    // 2) Get ISS position in ECI (km)
    Vec3 r_eci = issOrbitECI(t_sec);

    // 3) Rotate to ECEF (km)
    Vec3 r_ecef = eciToEcef(r_eci, t_sec);

    // 4) Convert to geodetic lat/lon/alt (all in km/degrees)
    double lon = std::atan2(r_ecef.y, r_ecef.x);
    double p   = std::hypot(r_ecef.x, r_ecef.y);

    constexpr double a = 6378.137, e2 = 6.69437999014e-3; // WGS84 (km)
    double θ    = std::atan2(r_ecef.z*a, p*(1-e2)*a);
    double sθ   = std::sin(θ), cθ = std::cos(θ);
    double lat  = std::atan2(
        r_ecef.z + e2*(1-e2)*a*sθ*sθ*sθ,
        p       - e2      *a*cθ*cθ*cθ
    );
    double N    = a / std::sqrt(1 - e2*std::sin(lat)*std::sin(lat));
    double alt  = p/std::cos(lat) - N;  // km

    // 5) Query IGRF in NED (north,east,down) in nT
    Vec3 B_ned = igrf(year, lat, lon, alt, "geodetic");

    // 6) Rotate NED➔ECEF➔ECI (all in nT)
    Vec3 B_ecef_nT = ned2eci(B_ned, lat * 180.0/M_PI, lon * 180.0/M_PI, t_sec);
    Vec3 B_eci_nT  = ecefToEci(B_ecef_nT, t_sec);

    // 7) Finally rotate ECI➔body, then convert nT➔T once
    Vec3 B_body_nT = rotateEciToBody(B_eci_nT, q);
    constexpr double nT2T = 1e-9;
    return B_body_nT * nT2T;
}
