// src/ned2eci.cpp
#include "ned2eci.h"
#include "Vec3.h"
#include <cmath>

// Earth’s sidereal rotation rate [rad/s]
static constexpr double omegaEarth = 7.2921150e-5;
static constexpr double deg2rad    = M_PI/180.0;
static constexpr double rad2deg    = 180.0/M_PI;

Vec3 ned2eci(const Vec3& ned,
             double latDeg,
             double lonDeg,
             double t_sec)
{
    // 1) NED → ECEF (still nT)
    double lat = latDeg * deg2rad;
    double lon = lonDeg * deg2rad;

    double Bn = ned.x;    // north
    double Be = ned.y;    // east
    double Bd = ned.z;    // down

    // rotation: first latitude about ECEF Y, then longitude about Z
    // equivalent to: [Bx_ecef;By_ecef;Bz_ecef] = R_z(lon) * R_y(lat) * [N;-D;E]
    double cosLat = std::cos(lat), sinLat = std::sin(lat);
    double cosLon = std::cos(lon), sinLon = std::sin(lon);

    // Note: down points toward Earth, but ECEF-z is up, so NED‐to‐ECEF is:
    //   X_e =  cosLon*( cosLat*(-Bd) ) - sinLon*Be - sinLat*(-Bd)*?
    // Simpler: rotate the vector [-Bd; Bn; Be] by latitude then longitude.
    // But here’s direct:
    double Bx_ecef =  cosLon*( cosLat*(-Bd) )
                     - sinLon*Be
                     - sinLat*Bn;
    double By_ecef =  sinLon*( cosLat*(-Bd) )
                     + cosLon*Be
                     - sinLat*Bn * 0; // actually north projection only on X/Z
    // Actually correct expansion from MATLAB version:
    Bx_ecef = cosLon * ( cosLat * -Bd ) - sinLat * Bn     - sinLon * Be;
    By_ecef = sinLon * ( cosLat * -Bd ) - sinLat * Bn     + cosLon * Be;
    double Bz_ecef = sinLat * ( -Bd ) + cosLat * Bn;

    // 2) ECEF → ECI rotation about Z by ω·t_sec
    double θ = omegaEarth * t_sec;
    double c = std::cos(θ), s = std::sin(θ);

    Vec3 ecef{ Bx_ecef, By_ecef, Bz_ecef };
    Vec3 eci {
        + c*ecef.x - s*ecef.y,
        + s*ecef.x + c*ecef.y,
          ecef.z
    };

    return eci;
}
