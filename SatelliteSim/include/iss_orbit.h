// include/iss_orbit.h
#pragma once

#include "Vec3.h"

/// Computes the satellite’s Earth-Centered Inertial (ECI) position
/// for a circular orbit at time t.
///
/// @param t  Time since start of simulation [s]
/// @returns  ECI position vector [x, y, z] in kilometers
Vec3 issOrbitECI(double t);
