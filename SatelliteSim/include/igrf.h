// include/igrf.h
#pragma once

#include "Vec3.h"
#include <vector>
#include <string>

/// Load and interpolate the IGRF g/h coefficient vector for a given decimal year.
/// Defined in igrf_loader.cpp
std::vector<double> loadIGRFCoeffs(double decimalYear);

/// Compute Earth’s magnetic field in the NED frame (north, east, down), in Tesla.
/// @param time    Decimal year (e.g. 2025.5)
/// @param lat     Geodetic latitude [rad]
/// @param lon     Geodetic longitude [rad]
/// @param alt     Altitude above WGS84 ellipsoid [m]
/// @param coord   “geodetic” (default) or “geocentric” coordinate mode
/// @returns       Magnetic field vector [Bn, Be, Bd] in Tesla
Vec3 igrf(double              time,
          double              lat,
          double              lon,
          double              alt,
          const std::string&  coord = "geodetic");
