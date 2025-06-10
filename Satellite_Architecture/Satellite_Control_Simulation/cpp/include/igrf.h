#ifndef IGRF_H
#define IGRF_H

#include <vector>
#include <string>

// Compute Earth's magnetic field via IGRF model.

// Scalar version:
//   - timeYear: decimal year (e.g. 2020.5)
//   - latitude_deg, longitude_deg: degrees
//   - altitude_km: if geodetic mode, height above surface; if geocentric, radius from center
//   - geodetic: true to treat (lat,alt) as geodetic, false if already geocentric
//   - Bx_nT, By_nT, Bz_nT: outputs in nanotesla
void igrf(double timeYear,
          double latitude_deg,
          double longitude_deg,
          double altitude_km,
          bool geodetic,
          double &Bx_nT,
          double &By_nT,
          double &Bz_nT);

// Vector version (reads .dat files from datDirectory):
//   - datDirectory: path where IGRF \"*.dat\" coefficient files reside
//   - timeYear: decimal year (e.g. 2020.5)
//   - latitude_deg, longitude_deg, altitude_km: vectors of same length
//   - geodetic: true if inputs are geodetic, false if geocentric
//   - Bx_nT, By_nT, Bz_nT: output vectors (same length), each component in nT
void igrf(const std::string &datDirectory,
          double timeYear,
          const std::vector<double> &latitude_deg,
          const std::vector<double> &longitude_deg,
          const std::vector<double> &altitude_km,
          bool geodetic,
          std::vector<double> &Bx_nT,
          std::vector<double> &By_nT,
          std::vector<double> &Bz_nT);

#endif // IGRF_H
