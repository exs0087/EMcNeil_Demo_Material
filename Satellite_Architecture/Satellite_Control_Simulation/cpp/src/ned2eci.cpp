#include "ned2eci.h"
#include <cmath>

std::array<double, 3> ned2eci(
    const std::array<double, 3> &ned_vec,
    double latitude_deg,
    double longitude_deg,
    double time_sec)
{
    // Convert degrees to radians
    constexpr double deg2rad = M_PI / 180.0;
    double lat_rad = latitude_deg * deg2rad;
    double lon_rad = longitude_deg * deg2rad;

    // Unpack the NED components (in nT)
    double BxNorth = ned_vec[0];
    double ByEast  = ned_vec[1];
    double BzDown  = ned_vec[2];

    // ----- NED → ECEF -----
    // Bx_ECEF = cos(lon) * [cos(lat)*(-BzDown)] - sin(lat)*BxNorth - sin(lon)*ByEast
    // By_ECEF = sin(lon) * [cos(lat)*(-BzDown)] - sin(lat)*BxNorth + cos(lon)*ByEast
    // Bz_ECEF = sin(lat) * (-BzDown) + cos(lat)*BxNorth

    double cos_lat = std::cos(lat_rad);
    double sin_lat = std::sin(lat_rad);
    double cos_lon = std::cos(lon_rad);
    double sin_lon = std::sin(lon_rad);

    double Bx_ECEF = cos_lon * (cos_lat * (-BzDown))
                   - sin_lat * BxNorth
                   - sin_lon * ByEast;
    double By_ECEF = sin_lon * (cos_lat * (-BzDown))
                   - sin_lat * BxNorth
                   + cos_lon * ByEast;
    double Bz_ECEF = sin_lat * (-BzDown)
                   + cos_lat * BxNorth;

    // ----- ECEF → ECI -----
    // Earth’s sidereal rotation: ω = 2π / 86164 rad/s
    constexpr double omega_dot_rad = 2.0 * M_PI / 86164.0; // [rad/s]
    double omega = omega_dot_rad * time_sec;              // [rad]

    double cos_omega = std::cos(omega);
    double sin_omega = std::sin(omega);

    // Rotate about the z‐axis:
    //   [Bx_ECI]   [ cosω  sinω   0 ] [Bx_ECEF]
    //   [By_ECI] = [-sinω  cosω   0 ] [By_ECEF]
    //   [Bz_ECI]   [   0     0    1 ] [Bz_ECEF]
    double Bx_ECI =  cos_omega * Bx_ECEF + sin_omega * By_ECEF;
    double By_ECI = -sin_omega * Bx_ECEF + cos_omega * By_ECEF;
    double Bz_ECI =  Bz_ECEF; // unchanged by yaw

    return {Bx_ECI, By_ECI, Bz_ECI};
}
