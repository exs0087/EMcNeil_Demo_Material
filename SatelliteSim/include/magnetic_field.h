// include/magnetic_field.h
#pragma once

#include "Vec3.h"
#include "Vec4.h"

/// Compute the magnetic field vector in the satellite body frame.
/// Internally:
///   1) Computes the ISS ECI position at time t
///   2) Converts to geodetic lat/lon/alt via ECEF
///   3) Queries IGRF for NED‐frame field
///   4) Transforms NED → ECEF → ECI
///   5) Rotates ECI vector into the body frame using the quaternion q
///
/// @param t  Time since simulation start [s]
/// @param q  Attitude quaternion (body→ECI, scalar-last)
/// @returns  Magnetic field in body frame [Bx, By, Bz] in Tesla
Vec3 getBodyFieldAt(double t, const Vec4& q);
