#ifndef MAGNETIC_FIELD_H
#define MAGNETIC_FIELD_H

#include "Vec3.h"
#include "Vec4.h"

/// Compute the Earth's magnetic field in the satellite body frame.
/// @param time    seconds since simulation start (for IGRF epoch lookup)
/// @param geo     geodetic coordinates:
///                geo.x = latitude  [rad]
///                geo.y = longitude [rad]
///                geo.z = altitude  [km]
///                geo.w = unused
/// @return        magnetic field in ECI frame, in Tesla
Vec3 getBodyFieldAt(double time, Vec4 const& geo);

#endif // MAGNETIC_FIELD_H
