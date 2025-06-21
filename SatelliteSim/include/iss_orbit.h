#pragma once

#include "Vec3.h"

// Earth’s rotation rate [rad/s]
inline constexpr double we = 7.2921150e-5;

// Propagates a circular, inclined orbit into ECEF
Vec3 issOrbitECI(double t);

void issOrbitGeod(double t, double& lat, double& lon, double& alt_km);
