#ifndef ISS_ORBIT_H
#define ISS_ORBIT_H

// Compute latitude, longitude, and altitude of a circular ISS-like orbit
// at time t (seconds). Assumes: inclination 52°, altitude 400 km,
// circular orbit (mu=398600 km^3/s^2, Re=6371 km).
//
// Inputs:
//   - t: elapsed time in seconds
// Outputs (by reference):
//   - latitude  [rad]  (geocentric; geodetic=geocentric since e=0)
//   - longitude [rad]  (in [0, 2*pi))
//   - altitude  [km]   (constant = 400)
void issOrbit(double t, double &latitude, double &longitude, double &altitude);

#endif // ISS_ORBIT_H
