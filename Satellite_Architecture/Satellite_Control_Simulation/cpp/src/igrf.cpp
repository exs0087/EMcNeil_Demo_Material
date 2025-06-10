// igrf.cpp
#include "igrf.h"
#include "igrf_loader.h"  // IGRFLoader
#include "legendre.hpp"    // computeLegendre
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <filesystem>

static constexpr double DEG2RAD = M_PI / 180.0;
static constexpr double Rearth_km = 6371.2;

// Helper: convert geodetic lat->geocentric colat and compute r, cd, sd
// Input: latitude_deg, altitude_km, geodetic flag
// Output: sinθ, cosθ, r, cd_multiplier, sd_multiplier
static void geodeticToGeocentric(double latitude_deg, double altitude_km,
                                 bool geodetic,
                                 double &sin_theta, double &cos_theta,
                                 double &r_km,
                                 double &cd, double &sd) {
    double lat_rad = latitude_deg * DEG2RAD;
    // Initial geocentric colatitude components:
    double costheta = std::cos((90.0 - latitude_deg) * DEG2RAD);
    double sintheta = std::sin((90.0 - latitude_deg) * DEG2RAD);
    if (geodetic) {
        // Convert from geodetic to geocentric
        double a = 6378.137;
        double f = 1.0 / 298.257223563;
        double b = a * (1.0 - f);
        double rho = std::hypot(a * sintheta, b * costheta);
        r_km = std::sqrt(altitude_km * altitude_km
                         + 2.0 * altitude_km * rho
                         + (a * a * a * a * sintheta * sintheta
                            + b * b * b * b * costheta * costheta) / (rho * rho));
        cd = (altitude_km + rho) / r_km;
        sd = ((a * a - b * b) / rho) * costheta * sintheta / r_km;
        double oldcost = costheta;
        // New geocentric cosθ, sinθ
        costheta = costheta * cd - sintheta * sd;
        sintheta = sintheta * cd + oldcost * sd;
    } else {
        // Geocentric provided directly: altitude_km is radius
        r_km = altitude_km;
        cd = 1.0;
        sd = 0.0;
    }
    sin_theta = sintheta;
    cos_theta = costheta;
}

// Core spherical-harmonic summation (general for multiple points)
static void computeField(const std::vector<double> &sin_theta,
                         const std::vector<double> &cos_theta,
                         const std::vector<double> &r_km,
                         const std::vector<double> &cd_vec,
                         const std::vector<double> &sd_vec,
                         const std::vector<double> &phi_rad,
                         bool geodetic,
                         IGRFLoader &loader,
                         double timeYear,
                         std::vector<double> &Bx_nT,
                         std::vector<double> &By_nT,
                         std::vector<double> &Bz_nT) {
    size_t Npts = sin_theta.size();
    if (phi_rad.size() != Npts || r_km.size() != Npts
        || cd_vec.size() != Npts || sd_vec.size() != Npts) {
        throw std::invalid_argument("Input vector sizes mismatch");
    }
    // Load interpolated g/h matrices
    int N = loader.maxDegree();
    // g, h are (N+1)x(N+1)
    std::vector<std::vector<double>> g(N+1, std::vector<double>(N+1, 0.0));
    std::vector<std::vector<double>> h(N+1, std::vector<double>(N+1, 0.0));
    loader.loadCoeffs(timeYear, g, h);

    // Pre-allocate field components and initialize to zero
    Bx_nT.assign(Npts, 0.0);
    By_nT.assign(Npts, 0.0);
    Bz_nT.assign(Npts, 0.0);
    std::vector<double> Br(Npts, 0.0);
    std::vector<double> Bt(Npts, 0.0);
    std::vector<double> Bp(Npts, 0.0);

    // Precompute cos(m*phi), sin(m*phi) for each point and m up to N
    std::vector<std::vector<double>> cos_m_phi(Npts, std::vector<double>(N+1));
    std::vector<std::vector<double>> sin_m_phi(Npts, std::vector<double>(N+1));
    for (size_t i = 0; i < Npts; ++i) {
        for (int m = 0; m <= N; ++m) {
            cos_m_phi[i][m] = std::cos(m * phi_rad[i]);
            sin_m_phi[i][m] = std::sin(m * phi_rad[i]);
        }
    }

    // For each point, precompute Legendre P[n][m] and dPdx[n][m] for all n,m
    // Store per-point vectors of 2D arrays
    std::vector<std::vector<std::vector<double>>> P_arr(Npts);
    std::vector<std::vector<std::vector<double>>> dPdx_arr(Npts);
    for (size_t i = 0; i < Npts; ++i) {
        P_arr[i].resize(N+1, std::vector<double>(N+1, 0.0));
        dPdx_arr[i].resize(N+1, std::vector<double>(N+1, 0.0));
        computeLegendre(sin_theta[i], P_arr[i], dPdx_arr[i]);
    }

    // Loop over degrees n and orders m
    for (int n = 1; n <= N; ++n) {
        // Precompute (Rearth/r)^(n+2) and (Rearth/r)^(n+1) per point
        std::vector<double> ar_n2(Npts), ar_n1(Npts);
        for (size_t i = 0; i < Npts; ++i) {
            double a_over_r = Rearth_km / r_km[i];
            ar_n2[i] = std::pow(a_over_r, n+2);
            ar_n1[i] = std::pow(a_over_r, n+1);
        }
        for (int m = 0; m <= n; ++m) {
            double g_nm = g[n][m];
            double h_nm = h[n][m];
            for (size_t i = 0; i < Npts; ++i) {
                double Pnm = P_arr[i][n][m];
                // dP/dθ = -cosθ * dP/dx
                double dP_dtheta = -cos_theta[i] * dPdx_arr[i][n][m];
                double term_cos = g_nm * cos_m_phi[i][m];
                double term_sin = h_nm * sin_m_phi[i][m];
                double common = (term_cos + term_sin);
                // Radial component
                Br[i] += (n + 1) * ar_n2[i] * common * Pnm;
                // Theta component (negative sign included below)
                Bt[i] -= ar_n1[i] * common * dP_dtheta;
                // Phi component
                if (m == 0) {
                    // Bp no contribution when m=0
                } else {
                    if (std::abs(sin_theta[i]) > 1e-16) {
                        Bp[i] -= ar_n1[i] * (m / sin_theta[i])
                                  * ( -g_nm * sin_m_phi[i][m] + h_nm * cos_m_phi[i][m] )
                                  * Pnm;
                    } else {
                        // At poles, use limit: 1/sinθ -> cotθ* derivative
                        Bp[i] -= ar_n1[i] * ( -g_nm * sin_m_phi[i][m] + h_nm * cos_m_phi[i][m] )
                                  * dP_dtheta;
                    }
                }
            }
        }
    }

    // Convert spherical (Br, Bt, Bp) -> Cartesian NED (Bx, By, Bz) per point
    for (size_t i = 0; i < Npts; ++i) {
        double bx = -Bt[i];
        double by = Bp[i];
        double bz = -Br[i];
        if (geodetic) {
            // Convert back to geodetic: Bx' = Bx*cd + Bz*sd; Bz' = Bz*cd - Bx*sd
            double bx_old = bx;
            bx = bx * cd_vec[i] + bz * sd_vec[i];
            bz = bz * cd_vec[i] - bx_old * sd_vec[i];
        }
        Bx_nT[i] = bx;
        By_nT[i] = by;
        Bz_nT[i] = bz;
    }
}

void igrf(double timeYear,
          double latitude_deg,
          double longitude_deg,
          double altitude_km,
          bool geodetic,
          double &Bx_nT,
          double &By_nT,
          double &Bz_nT) {
    // Scalar wrapper: convert scalars to single-element vectors
    std::vector<double> lat{latitude_deg}, lon{longitude_deg}, alt{altitude_km};
    std::vector<double> Bx, By, Bz;
    // Use current directory for dat files
    std::string datDir = "./datfiles/igrf_all_epochs.json";
    igrf(datDir, timeYear, lat, lon, alt, geodetic, Bx, By, Bz);
    Bx_nT = Bx[0];
    By_nT = By[0];
    Bz_nT = Bz[0];
}

void igrf(const std::string &datDirectory,
          double timeYear,
          const std::vector<double> &latitude_deg,
          const std::vector<double> &longitude_deg,
          const std::vector<double> &altitude_km,
          bool geodetic,
          std::vector<double> &Bx_nT,
          std::vector<double> &By_nT,
          std::vector<double> &Bz_nT) {
    size_t Np = latitude_deg.size();
    if (longitude_deg.size() != Np || altitude_km.size() != Np) {
        throw std::invalid_argument("Latitude, longitude, and altitude must have same size");
    }
    // Precompute per-point quantities
    std::vector<double> sin_theta(Np), cos_theta(Np), r_km(Np), cd_vec(Np), sd_vec(Np), phi_rad(Np);
    for (size_t i = 0; i < Np; ++i) {
        geodeticToGeocentric(latitude_deg[i], altitude_km[i], geodetic,
                              sin_theta[i], cos_theta[i], r_km[i], cd_vec[i], sd_vec[i]);
        phi_rad[i] = longitude_deg[i] * DEG2RAD;
    }
    // Initialize loader and compute field
    IGRFLoader loader(datDirectory);
    computeField(sin_theta, cos_theta, r_km, cd_vec, sd_vec, phi_rad, geodetic,
                 loader, timeYear, Bx_nT, By_nT, Bz_nT);
}
