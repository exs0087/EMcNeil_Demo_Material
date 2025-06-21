// src/igrf.cpp

#include "igrf.h"            // declares Vec3, geodeticToGeocentric(), computeSphericalField(), igrf()
#include "igrf_loader.h"     // loadIGRFCoeffs()
#include <cmath>
#include <tuple>

// Reference radius for IGRF (km)
static constexpr double Re_IGRF = 6371.2;

// WGS-84 ellipsoid constants (km)
static constexpr double a_ell = 6378.137;
static constexpr double f_ell = 1.0/298.257223563;
static constexpr double b_ell = a_ell*(1 - f_ell);

//------------------------------------------------------------------------------
// Convert geodetic (lat,lon,alt) → geocentric spherical (r,θ,φ)
//   lat,lon in radians; alt in km
// Outputs:
//   r_km      – radius from Earth center [km]
//   cos_theta – cos(colatitude) = z/r
//   sin_theta – sin(colatitude)
//------------------------------------------------------------------------------
void geodeticToGeocentric(double lat_rad,
                          double lon_rad,
                          double alt_km,
                          double &r_km,
                          double &cos_theta,
                          double &sin_theta)
{
    double lat_deg   = lat_rad * 180.0/M_PI;
    double costheta0 = std::cos((90.0 - lat_deg) * M_PI/180.0);
    double sintheta0 = std::sin((90.0 - lat_deg) * M_PI/180.0);

    double rho = std::hypot(a_ell * sintheta0,
                            b_ell * costheta0);

    r_km = std::sqrt(
        alt_km*alt_km
      + 2.0*alt_km*rho
      + ((a_ell*a_ell*a_ell*a_ell)*sintheta0*sintheta0
       + (b_ell*b_ell*b_ell*b_ell)*costheta0*costheta0) / (rho*rho)
    );

    double cd = (alt_km + rho) / r_km;
    double sd = ((a_ell*a_ell - b_ell*b_ell) / rho) * (costheta0*sintheta0) / r_km;

    cos_theta = costheta0*cd - sintheta0*sd;
    sin_theta = sintheta0*cd + costheta0*sd;
}

//------------------------------------------------------------------------------
// Compute the spherical-harmonic magnetic field at radius r_km, colatitude
// theta_rad, longitude phi_rad, using the flat GH coefficient vector.
// Returns (Br [nT], Bt [nT], Bp [nT]) in the spherical basis.
//------------------------------------------------------------------------------
std::tuple<double,double,double>
computeSphericalField(double r_km,
                      double theta_rad,
                      double phi_rad,
                      const std::vector<double> &gh)
{
    // derive maximum degree nmax
    int nmax = static_cast<int>(std::sqrt(gh.size() + 1) - 1);

    // precompute sin/cos of colatitude & longitude
    double ct = std::cos(theta_rad), st = std::sin(theta_rad);
    double pr = Re_IGRF / r_km;

    // storage for P(n,m) and dP(n,m)
    std::vector<std::vector<double>> P(nmax+1, std::vector<double>(nmax+1, 0.0));
    std::vector<std::vector<double>> dP(nmax+1, std::vector<double>(nmax+1, 0.0));

    // initial
    P[0][0] = 1.0;
    dP[0][0] = 0.0;

    // build associated Legendre
    for (int n = 1; n <= nmax; ++n) {
        for (int m = 0; m <= n; ++m) {
            if (n == m) {
                double k = std::sqrt(1.0 - 1.0/(2.0*n));
                P[n][m]   = k * st * P[n-1][m-1];
                dP[n][m]  = k * (ct * P[n-1][m-1] + st * dP[n-1][m-1]);
            }
            else if (n == 1 && m == 0) {
                P[n][m]   = ct * P[n-1][m];
                dP[n][m]  = -st * P[n-1][m] + ct * dP[n-1][m];
            }
            else {
                double a = ((2.0*n - 1.0)/(n - m));
                double b = ((n + m - 1.0)/(n - m));
                P[n][m]   = a * ct * P[n-1][m] - b * P[n-2][m];
                dP[n][m]  = a * (ct * dP[n-1][m] - st * P[n-1][m])
                          - b * dP[n-2][m];
            }
        }
    }

    double Br = 0.0, Bt = 0.0, Bp = 0.0;
    int idx = 0;
    for (int n = 1; n <= nmax; ++n) {
        double ar = std::pow(pr, n+2);
        for (int m = 0; m <= n; ++m) {
            double gnm = gh[idx++];
            double hnm = (m>0 ? gh[idx++] : 0.0);
            double cosm = std::cos(m*phi_rad);
            double sinm = std::sin(m*phi_rad);
            double coeff = gnm * cosm + hnm * sinm;

            Br += (n+1) * ar * coeff * P[n][m];
            Bt -= ar * coeff * dP[n][m];
            if (m>0) {
                double coeff2 = -gnm * sinm + hnm * cosm;
                Bp += ar * (m * coeff2 * P[n][m] / st);
            }
        }
    }

    return {Br, Bt, Bp};
}

//------------------------------------------------------------------------------
// Main IGRF entry: geodetic→NED magnetic field [nT].
//------------------------------------------------------------------------------
Vec3 igrf(double time,
          double lat_rad,
          double lon_rad,
          double alt_km,
          const std::string & /*coord*/)
{
    double r_km, cos_theta, sin_theta;
    geodeticToGeocentric(lat_rad, lon_rad, alt_km,
                         r_km, cos_theta, sin_theta);

    double theta = std::atan2(sin_theta, cos_theta);
    double phi   = lon_rad;

    auto gh = loadIGRFCoeffs(time);
    auto [Br, Bt, Bp] = computeSphericalField(r_km, theta, phi, gh);

    return Vec3{ -Bt,   // north
                  Bp,   // east
                 -Br }; // down
}
