// src/magnetic_field.cpp

#include "magnetic_field.h"
#include "igrf.h"        // igrf()
#include <cmath>

// Earth's rotation rate (rad/s)
static constexpr double we = 7.2921150e-5;

// NED → ECEF rotation (leave as is)
static Vec3 ned2ecef(const Vec3& v_ned, double lat, double lon) {
    double s_lat = std::sin(lat), c_lat = std::cos(lat);
    double s_lon = std::sin(lon), c_lon = std::cos(lon);
    return Vec3{
        -s_lat * c_lon * v_ned.x + -s_lat * s_lon * v_ned.y +  c_lat * v_ned.z,
        -       s_lon * v_ned.x +      c_lon * v_ned.y +  0.0       * v_ned.z,
        -c_lat * c_lon * v_ned.x + -c_lat * s_lon * v_ned.y + -s_lat * v_ned.z
    };
}

Vec3 getBodyFieldAt(double time, Vec4 const& geo) {
    // geo.x = geodetic latitude [rad]
    // geo.y = longitude           [rad]
    // geo.z = altitude above ellipsoid [km]

    // 1) IGRF returns NED in nT (nanotesla)
    Vec3 B_ned_nT = igrf(time, geo.x, geo.y, geo.z, "geodetic");

    // 2) Convert nT → T for use everywhere downstream
    Vec3 B_ned{ B_ned_nT.x * 1e-9,
                B_ned_nT.y * 1e-9,
                B_ned_nT.z * 1e-9 };

    // 3) Rotate NED → ECEF (now in Tesla)
    Vec3 B_ecef = ned2ecef(B_ned, geo.x, geo.y);

    // 4) Rotate ECEF → ECI by undoing Earth rotation
    double gst = we * time;
    double c   = std::cos(gst);
    double s   = std::sin(gst);
    Vec3 B_eci{
        c * B_ecef.x + s * B_ecef.y,
       -s * B_ecef.x + c * B_ecef.y,
        B_ecef.z
    };

    // 5) Return in ECI frame, in **Tesla**
    return B_eci;
}
