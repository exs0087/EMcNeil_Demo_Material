// ode4m.cpp
#include "ode4m.h"
#include "eulerseqns2.h" // provides extern std::array<double,3> L and function eulerseqns2(...)
#include "iss_orbit.h"   // issOrbit
#include "igrf.h"        // igrf
#include "ned2eci.h"     // ned2eci

#include <cmath>
#include <stdexcept>
#include <array>
#include <vector>
#include <algorithm>

enum { NX = 7 }; // state dimension: [w1,w2,w3,q1,q2,q3,q4]

// Global history vectors for plotting/analysis:
std::vector<std::array<double,3>> LHIST;  // control torque history (body‐frame)
std::vector<std::array<double,3>> MHIST;  // dipole moment history (body‐frame)

std::vector<std::array<double,7>> ode4m(
    const std::function<std::array<double,7>(double, const std::array<double,7>&)> &odefun,
    const std::vector<double> &tspan,
    const std::array<double,7> &y0,
    int cmode
) {
    size_t N = tspan.size();
    if (N < 2) {
        throw std::invalid_argument("tspan must have at least two entries");
    }
    // Step sizes h[i] = tspan[i] - tspan[i-1]
    std::vector<double> h(N);
    for (size_t i = 1; i < N; ++i) {
        h[i] = tspan[i] - tspan[i-1];
        if (h[i] <= 0) {
            throw std::invalid_argument("tspan entries must be strictly increasing");
        }
    }
    // Pre-allocate solution array: Y[i] is state at tspan[i]
    std::vector<std::array<double,7>> Y(N);
    Y[0] = y0;

    // Initialize history arrays
    LHIST.assign(N, {0.0,0.0,0.0});
    MHIST.assign(N, {0.0,0.0,0.0});

    // Pre-simulate Earth's magnetic field at 1 Hz for 16001 seconds
    const size_t M = 16001;
    std::vector<double> Bx_i(M), By_i(M), Bz_i(M);
    for (size_t step = 0; step < M; ++step) {
        double dt = static_cast<double>(step);
        double lat, lon, alt;
        issOrbit(dt, lat, lon, alt); // lat,lon in rad, alt in km
        // IGRF wants lat/lon in degrees, altitude in km
        double Bxn, Byn, Bzn;
        igrf(dt,                   // time (we could pass real epoch here if needed)
             lat * 180.0/M_PI,     // geodetic latitude [deg]
             lon * 180.0/M_PI,     // geodetic longitude [deg]
             alt,                  // altitude [km]
             true,                 // geodetic flag
             Bxn, Byn, Bzn);       // output in nT, in NED
        std::array<double,3> ned = { Bxn, Byn, Bzn };
        auto eci = ned2eci(ned,
                          lat * 180.0/M_PI,
                          lon * 180.0/M_PI,
                          dt);
        Bx_i[step] = eci[0];
        By_i[step] = eci[1];
        Bz_i[step] = eci[2];
    }

    // Spacecraft inertia properties (same rectangular prism)
    double mass = 10.0; // kg
    double height = 0.34, width = 0.20, depth = 0.10; // m
    double Jx = mass/12.0 * (width*width + depth*depth);
    double Jy = mass/12.0 * (height*height + depth*depth);
    double Jz = mass/12.0 * (height*height + width*width);
    double d1 = 1.0 / Jx;
    double d2 = 1.0 / Jy;
    double d3 = 1.0 / Jz;

    // Control parameters
    const double moment = 0.2;     // A·m² (max dipole)
    const double sample_rate = 1.0; // sec

    // Variables for B-dot or ω×B logic
    double tlast = 0.0;
    double Bxlast = 0.0, Bylast = 0.0, Bzlast = 0.0;
    double mxlast = 0.0, mylast = 0.0, mzlast = 0.0;

    // Temporary arrays for RK4 and rotation
    std::array<double,7> F1, F2, F3, F4;
    std::array<double,7> yi, ytemp;

    for (size_t i = 1; i < N; ++i) {
        double ti = tspan[i-1];
        double hi = h[i];
        yi = Y[i-1];

        // 1) Extract angular rates & quaternion from yi
        double w1 = yi[0], w2 = yi[1], w3 = yi[2];
        double q1 = yi[3], q2 = yi[4], q3 = yi[5], q4 = yi[6];

        // Build body←ECI DCM from quaternion:
        double C11 = 1 - 2*(q2*q2 + q3*q3);
        double C12 = 2*(q1*q2 + q3*q4);
        double C13 = 2*(q3*q1 - q2*q4);
        double C21 = 2*(q2*q1 - q3*q4);
        double C22 = 1 - 2*(q3*q3 + q1*q1);
        double C23 = 2*(q2*q3 + q1*q4);
        double C31 = 2*(q3*q1 + q2*q4);
        double C32 = 2*(q2*q3 - q1*q4);
        double C33 = 1 - 2*(q1*q1 + q2*q2);

        // 2) Interpolate Earth‐field at time ti (linear interp if ti not integer)
        double idx_f = std::floor(ti);
        size_t i0 = (idx_f < 0) ? 0 : static_cast<size_t>(idx_f);
        if (i0 >= M-1) i0 = M-1;
        size_t i1 = std::min(i0 + 1, M-1);
        double frac = ti - std::floor(ti);
        double Bxeci = (1 - frac)*Bx_i[i0] + frac*Bx_i[i1];
        double Byeci = (1 - frac)*By_i[i0] + frac*By_i[i1];
        double Bzeci = (1 - frac)*Bz_i[i0] + frac*Bz_i[i1];

        // 3) Rotate ECI‐field into body frame: B_body = DCM_body←ECI * B_ECI
        double Bx_body = C11*Bxeci + C12*Byeci + C13*Bzeci;
        double By_body = C21*Bxeci + C22*Byeci + C23*Bzeci;
        double Bz_body = C31*Bxeci + C32*Byeci + C33*Bzeci;
        // Convert from nT → Tesla
        std::array<double,3> B_body_vec = { Bx_body*1e-9,
                                            By_body*1e-9,
                                            Bz_body*1e-9 };

        // 4) Compute control dipole (A·m²) based on cmode
        std::array<double,3> mvec_body = {0.0, 0.0, 0.0};

        if (cmode == 1) {
            // ---- B-dot (bang-bang) control ----
            if (std::fmod(ti, sample_rate) < 1e-6) {
                if (std::fabs(ti - tlast) > 1e-8) {
                    double Bdotx = (Bx_body - Bxlast) / (ti - tlast);
                    double Bdoty = (By_body - Bylast) / (ti - tlast);
                    double Bdotz = (Bz_body - Bzlast) / (ti - tlast);
                    mvec_body = {
                        -moment * ((Bdotx > 0) - (Bdotx < 0)),
                        -moment * ((Bdoty > 0) - (Bdoty < 0)),
                        -moment * ((Bdotz > 0) - (Bdotz < 0))
                    };
                    Bxlast = Bx_body;
                    Bylast = By_body;
                    Bzlast = Bz_body;
                    tlast = ti;
                } else {
                    // First sample at t=0
                    mvec_body = {0.0, 0.0, 0.0};
                    Bxlast = Bx_body;
                    Bylast = By_body;
                    Bzlast = Bz_body;
                    tlast = ti;
                }
            } else {
                // Hold previous m
                mvec_body = { mxlast, mylast, mzlast };
            }
        }
        else if (cmode == 2) {
            // ---- ω × B “damping” control ----
            if (std::fmod(ti, sample_rate) < 1e-6) {
                if (std::fabs(ti - tlast) > 1e-8) {
                    double Bnorm = std::sqrt(
                        B_body_vec[0]*B_body_vec[0] +
                        B_body_vec[1]*B_body_vec[1] +
                        B_body_vec[2]*B_body_vec[2]
                    );
                    // Normalized field unit vector
                    std::array<double,3> b_unit = {
                        B_body_vec[0]/Bnorm,
                        B_body_vec[1]/Bnorm,
                        B_body_vec[2]/Bnorm
                    };
                    // Cross product w × b
                    std::array<double,3> w_vec = { w1, w2, w3 };
                    std::array<double,3> cross_wb = {
                        w_vec[1]*b_unit[2] - w_vec[2]*b_unit[1],
                        w_vec[2]*b_unit[0] - w_vec[0]*b_unit[2],
                        w_vec[0]*b_unit[1] - w_vec[1]*b_unit[0]
                    };
                    double k = 1.0; // scaling factor
                    mvec_body = {
                        (k / Bnorm) * cross_wb[0],
                        (k / Bnorm) * cross_wb[1],
                        (k / Bnorm) * cross_wb[2]
                    };
                    // Clip each component to ±0.2 A·m²
                    for (int c=0; c<3; ++c) {
                        if (mvec_body[c] > moment)       mvec_body[c] = moment;
                        else if (mvec_body[c] < -moment) mvec_body[c] = -moment;
                    }
                    Bxlast = Bx_body;
                    Bylast = By_body;
                    Bzlast = Bz_body;
                    tlast = ti;
                } else {
                    // First sample at t=0
                    mvec_body = {0.0, 0.0, 0.0};
                    Bxlast = Bx_body;
                    Bylast = By_body;
                    Bzlast = Bz_body;
                    tlast = ti;
                }
            } else {
                // Hold previous m
                mvec_body = { mxlast, mylast, mzlast };
            }
        }

        // Remember this dipole for the next iteration
        mxlast = mvec_body[0];
        mylast = mvec_body[1];
        mzlast = mvec_body[2];

        // 5) Compute control torque in body frame: τ = m × B
        std::array<double,3> torque_body = {
            mvec_body[1]*B_body_vec[2] - mvec_body[2]*B_body_vec[1],
            mvec_body[2]*B_body_vec[0] - mvec_body[0]*B_body_vec[2],
            mvec_body[0]*B_body_vec[1] - mvec_body[1]*B_body_vec[0]
        };
        // Scale by 1/J to feed into Euler’s eqn
        L[0] = torque_body[0] * d1;
        L[1] = torque_body[1] * d2;
        L[2] = torque_body[2] * d3;

        // Store into history arrays
        LHIST[i] = { L[0], L[1], L[2] };
        MHIST[i] = mvec_body;

        // 6) Classical RK4 (call eulerseqns2 via the lambda 'odefun'):
        F1 = odefun(ti, yi);

        for (int j = 0; j < NX; ++j) {
            ytemp[j] = yi[j] + 0.5 * hi * F1[j];
        }
        F2 = odefun(ti + 0.5*hi, ytemp);

        for (int j = 0; j < NX; ++j) {
            ytemp[j] = yi[j] + 0.5 * hi * F2[j];
        }
        F3 = odefun(ti + 0.5*hi, ytemp);

        for (int j = 0; j < NX; ++j) {
            ytemp[j] = yi[j] + hi * F3[j];
        }
        F4 = odefun(tspan[i], ytemp);

        for (int j = 0; j < NX; ++j) {
            Y[i][j] = yi[j]
                    + (hi / 6.0) * (F1[j] + 2.0*F2[j] + 2.0*F3[j] + F4[j]);
        }

        // 7) **Normalize quaternion** so ‖q‖ = 1
        double rq1 = Y[i][3], rq2 = Y[i][4], rq3 = Y[i][5], rq4 = Y[i][6];
        double normq = std::sqrt(rq1*rq1 + rq2*rq2 + rq3*rq3 + rq4*rq4);
        Y[i][3] = rq1 / normq;
        Y[i][4] = rq2 / normq;
        Y[i][5] = rq3 / normq;
        Y[i][6] = rq4 / normq;
    }

    return Y;
}
