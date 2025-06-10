#include "control_sim.hpp"
#include <cmath>
#include <algorithm>

// element-wise scale and add for our 7-vector state (w,q)
static Vec7 vec7_scale(const Vec7 &v, double s) {
    Vec7 out;
    for (int i = 0; i < 7; ++i) out[i] = v[i] * s;
    return out;
}
static Vec7 vec7_add(const Vec7 &a, const Vec7 &b) {
    Vec7 out;
    for (int i = 0; i < 7; ++i) out[i] = a[i] + b[i];
    return out;
}

// simple 3×3 cross‐product
static Vec3 cross3(const Vec3 &a, const Vec3 &b) {
    return Vec3{
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

History simulateControl(double t0,
                        double tf,
                        double dt,
                        const Vec3 &w0,
                        const Vec4 &q0,
                        int cmode) 
{
    // build time‐vector
    size_t N = size_t((tf - t0)/dt) + 1;
    std::vector<double> tspan(N);
    for (size_t i = 0; i < N; ++i) tspan[i] = t0 + i*dt;

    // initial state [ω; q]
    Vec7 y;
    y[0]=w0[0]; y[1]=w0[1]; y[2]=w0[2];
    y[3]=q0[0]; y[4]=q0[1]; y[5]=q0[2]; y[6]=q0[3];

    // storage
    std::vector<Vec7> Y(N);
    std::vector<Vec3> Lhist(N);
    std::vector<Vec3> Mhist(N);
    Y[0] = y;

    // constants from your MATLAB script
    const double moment = 0.2;  
    const double sample_rate = 1.0;
    // spacecraft inertias (matching euler_dynamics.cpp)
    constexpr double mass = 10, h = 0.34, w = 0.2, d = 0.1;
    constexpr double Jx = mass / 12 * (w * w + d * d),
                     Jy = mass / 12 * (h * h + d * d),
                     Jz = mass / 12 * (h * h + w * w);
    const double c1 = (Jy - Jz) / Jx;
    const double d1 = 1.0 / Jx;
    const double c2 = (Jz - Jx) / Jy;
    const double d2 = 1.0 / Jy;
    const double c3 = (Jx - Jy) / Jz;
    const double d3 = 1.0 / Jz;

    // previous‐sample storage
    double tlast = 0, Bxlast=0, Bylast=0, Bzlast=0;
    Vec3 mxlast{0,0,0};
    Vec3 mylast{0,0,0};
    Vec3 mzlast{0,0,0};

    // loop over steps
    for (size_t i = 1; i < N; ++i) {
        double ti = tspan[i-1];
        double hi = dt;
        Vec7 yi = Y[i-1];

        // unpack state
        Vec3 w{ yi[0], yi[1], yi[2] };
        Vec4 q{ yi[3], yi[4], yi[5], yi[6] };

        // ——— controller section ———
        Vec3 torque_body;
        Vec3 mvec_body;

        // first, compute body‐frame B from ECI (you’ll already have this callback)
        Vec3 Beci = getBodyFieldAt(ti);   // you must implement this
        Vec3 Bvec = Beci * 1e-9;          // nT→T

        if (cmode == 1) {
            // bang-bang B-dot control (exactly as in ode4m.m) :contentReference[oaicite:2]{index=2}:contentReference[oaicite:3]{index=3}
            double Bdotx = (Bvec[0] - Bxlast)/(ti - tlast);
            double Bdoty = (Bvec[1] - Bylast)/(ti - tlast);
            double Bdotz = (Bvec[2] - Bzlast)/(ti - tlast);

            mvec_body = Vec3{
                -moment * std::copysign(1.0, Bdotx),
                -moment * std::copysign(1.0, Bdoty),
                -moment * std::copysign(1.0, Bdotz)
            };
        }
        else if (cmode == 2) {
            // ω×B controller :contentReference[oaicite:4]{index=4}:contentReference[oaicite:5]{index=5}
            double Bnorm = std::sqrt(Bvec[0]*Bvec[0]
                                    + Bvec[1]*Bvec[1]
                                    + Bvec[2]*Bvec[2]);
            Vec3 bdir{ Bvec[0]/Bnorm,
                       Bvec[1]/Bnorm,
                       Bvec[2]/Bnorm };

            double k = 1.0;
            // raw dipole
            Vec3 raw{ w[1]*bdir[2] - w[2]*bdir[1],
                      w[2]*bdir[0] - w[0]*bdir[2],
                      w[0]*bdir[1] - w[1]*bdir[0] };
            for (auto &x : raw) x *= k / Bnorm;

            // saturate to ±moment
            mvec_body = Vec3{
                std::clamp(raw[0], -moment, moment),
                std::clamp(raw[1], -moment, moment),
                std::clamp(raw[2], -moment, moment)
            };
        }

        // common cross‐product to get torque
        torque_body = cross3(mvec_body, Bvec);

        // store histories
        Lhist[i] = Vec3{ torque_body[0]*d1,
                         torque_body[1]*d2,
                         torque_body[2]*d3 };
        Mhist[i] = mvec_body;

        // update last‐sample values
        tlast = ti;
        Bxlast = Bvec[0];
        Bylast = Bvec[1];
        Bzlast = Bvec[2];

        // ——— integrate one RK4 step ———
        Vec7 k1 = eulerDynamics(ti,           yi,           Lhist[i]);
        Vec7 k2 = eulerDynamics(ti+hi/2,      vec7_add(yi, vec7_scale(k1,hi/2)), Lhist[i]);
        Vec7 k3 = eulerDynamics(ti+hi/2,      vec7_add(yi, vec7_scale(k2,hi/2)), Lhist[i]);
        Vec7 k4 = eulerDynamics(tspan[i],     vec7_add(yi, vec7_scale(k3,hi)),   Lhist[i]);

        Vec7 ynew = yi;
        for (int j = 0; j < 7; ++j)
            ynew[j] += (hi/6.0)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

        Y[i] = ynew;
    }

    // wrap up
    History H;
    H.t = std::move(tspan);
    H.states = std::move(Y);
    H.torques = std::move(Lhist);
    H.mdipoles = std::move(Mhist);
    return H;
}
