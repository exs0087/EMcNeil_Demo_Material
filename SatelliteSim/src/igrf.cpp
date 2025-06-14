// src/igrf.cpp

#include "igrf.h"              // declares Vec3 igrf(...) and loadIGRFCoeffs(...)
#include <cmath>
#include <stdexcept>
#include <vector>              // for std::vector
#include <string>              // for std::string_view if you switch

static constexpr double pi      = 3.141592653589793;
static constexpr double deg2rad = pi/180.0;
static constexpr double rad2deg = 180.0/pi;
static constexpr double Re_km   = 6371.2;      // IGRF reference radius [km]

Vec3 igrf(double time,      // decimal year
          double latRad,    // geodetic latitude [rad]
          double lonRad,    // geodetic longitude [rad]
          double altM,      // altitude [m]
          const std::string& coord /*="geodetic"*/)
{
    // 1) Convert units
    double altKm = altM * 1e-3;
    double latDeg = latRad * rad2deg;
    double lonDeg = lonRad * rad2deg;

    // 2) Geodetic ↔ Geocentric conversion
    const double a = 6378.137, f = 1.0/298.257223563;
    const double b = a * (1 - f);

    double costheta = std::cos((90.0 - latDeg) * deg2rad);
    double sintheta = std::sin((90.0 - latDeg) * deg2rad);
    double cd = 1.0, sd = 0.0, r;

    if (coord == "geodetic") {
        double rho = std::hypot(a * sintheta, b * costheta);
        r = std::sqrt(altKm*altKm
                    + 2*altKm*rho
                    + (a*a*sintheta*sintheta + b*b*costheta*costheta)/(rho*rho));
        cd = (altKm + rho) / r;
        sd = ((a*a - b*b) / rho) * costheta * sintheta / r;
        // update costheta/sintheta
        double old = costheta;
        costheta = costheta*cd - sintheta*sd;
        sintheta = old*sd     + sintheta*cd;
    }
    else if (coord == "geocentric") {
        r = altKm;
    }
    else {
        throw std::invalid_argument("igrf: invalid coord flag");
    }

    // 3) Load & interpolate the g/h coefficients
    auto gh = loadIGRFCoeffs(time);
    int nmax = static_cast<int>(std::sqrt(gh.size() + 1) - 1);

    // 4) Precompute longitude terms
    std::vector<double> cosphi(nmax+1), sinphi(nmax+1);
    for (int m = 0; m <= nmax; ++m) {
        cosphi[m] = std::cos(m * lonRad);
        sinphi[m] = std::sin(m * lonRad);
    }

    // 5) Build Schmidt‐normalized Legendre & derivative arrays
    int Psize = (nmax+1)*(nmax+2)/2;
    std::vector<double> P(Psize), dP(Psize);
    P[0] = 1.0;
    dP[0] = 0.0;
    for (int n = 1; n <= nmax; ++n) {
        for (int m = 0; m <= n; ++m) {
            int idx = n*(n+1)/2 + m;
            if (n == m) {
                double fac = std::sqrt(1.0 - 1.0/(2.0*n));
                P[idx]  = fac * sintheta * P[idx - n - 1];
                dP[idx] = fac * (sintheta * dP[idx - n - 1]
                               + costheta * P[idx - n - 1]);
            } else {
                double k1 = (2.0*n - 1) / std::sqrt(n*n - m*m);
                double k2 = std::sqrt(((n-1.0)*(n-1.0) - m*m)/(n*n - m*m));
                P[idx]  = k1*(costheta*P[idx-n] - sintheta*P[idx-2*n-1])
                         - k2 * P[idx-2*n-1];
                dP[idx] = k1*(costheta*dP[idx-n] - sintheta*dP[idx-2*n-1])
                         - k2 * dP[idx-2*n-1];
            }
        }
    }

    // 6) Spherical‐harmonic summation
    double Br = 0.0, Bt = 0.0, Bp = 0.0;
    double ar = std::pow(Re_km / r, 2);
    int coeffIdx = 0;
    for (int n = 1; n <= nmax; ++n) {
        ar *= (Re_km / r);  // now (Re/r)^(n+2)
        for (int m = 0; m <= n; ++m) {
            double gnm = gh[coeffIdx++];
            double hnm = (m > 0 ? gh[coeffIdx++] : 0.0);
            int idx  = n*(n+1)/2 + m;
            double term = gnm*cosphi[m] + hnm*sinphi[m];

            Br += ar * (n+1) * term * P[idx];
            Bt -= ar * term * dP[idx];
            if (m > 0) {
                // handle the 1/sintheta vs dP case
                double fac = (sintheta != 0.0 ? P[idx] : dP[idx]) / (sintheta != 0.0 ? sintheta : 1.0);
                Bp -= ar * m * ( -gnm*sinphi[m] + hnm*cosphi[m] ) * fac;
            }
        }
    }

    // 7) Convert spherical->NED & apply geodetic correction
    double Bn = -Bt;
    double Be =  Bp;
    double Bd = -Br;

    double oldBn = Bn;
    Bn = Bn*cd + Bd*sd;
    Bd = Bd*cd - oldBn*sd;

    // 8) Convert nT -> T
    constexpr double nT2T = 1e-9;
    return Vec3{ Bn * nT2T, Be * nT2T, Bd * nT2T };
}
