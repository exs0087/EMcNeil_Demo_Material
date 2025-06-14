// include/ned2eci.h
#pragma once

#include "Vec3.h"

/// Transform a magnetic field vector from NED (north, east, down) frame to ECI frame.
/// Internally:
///   1) Rotates NED→ECEF based on geodetic lat/lon
///   2) Rotates ECEF→ECI by Earth’s sidereal rotation at time t
///
/// @param ned   Magnetic field vector [Bn, Be, Bd] in Tesla (NED frame)
/// @param lat   Geodetic latitude [rad]
/// @param lon   Geodetic longitude [rad]
/// @param time  Seconds since simulation start (for Earth rotation)
/// @returns     Magnetic field vector [x, y, z] in ECI frame (Tesla)
Vec3 ned2eci(const Vec3& ned, double lat, double lon, double time);
