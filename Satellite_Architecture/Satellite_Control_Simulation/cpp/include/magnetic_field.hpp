#pragma once
#include <array>
using Vec3 = std::array<double,3>;

/// Example stubs—port your MATLAB routines here:
Vec3 issOrbitECI(double t);
Vec3 igrfECEF(double t, double lat, double lon, double alt);
