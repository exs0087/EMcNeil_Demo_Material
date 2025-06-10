#ifndef NED2ECI_H
#define NED2ECI_H

#include <array>

// Convert a magnetic field vector from NED (North, East, Down) to ECI.
//
// Inputs:
//   ned_vec       – {BxNorth [nT], ByEast [nT], BzDown [nT]}
//   latitude_deg  – geodetic latitude in degrees
//   longitude_deg – geodetic longitude in degrees
//   time_sec      – seconds since the reference epoch (e.g., UT1=0 for Earth rotation)
//
// Output (returned array):
//   {Bx_ECI, By_ECI, Bz_ECI} in nanotesla
std::array<double, 3> ned2eci(
    const std::array<double, 3> &ned_vec,
    double latitude_deg,
    double longitude_deg,
    double time_sec);

#endif // NED2ECI_H
