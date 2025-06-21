// src/control_algorithm.cpp

#include <cmath>
#include <iostream>
#include "control_algorithm.h"   // for ControlMode, SensorReadings, ActuatorCommands

// Default to proportional ω×B mode (no auto‐switch)
static ControlMode g_mode = CTRL_MODE_WX_B;

extern "C" void control_algorithm_set_mode(ControlMode mode) {
    g_mode = mode;
}

extern "C" void control_algorithm(const SensorReadings* sens,
                                  ActuatorCommands*   act)
{
    // ——— DEBUG: See what branch you’re in, Bnorm and ω×b ———
    // Bx, By, Bz are **always in Tesla** here!
    if (sens->t < 10.0) {
        double Bnorm = std::sqrt(sens->Bx*sens->Bx
                               + sens->By*sens->By
                               + sens->Bz*sens->Bz);
        double cx =  sens->wy * sens->Bz - sens->wz * sens->By;
        double cy =  sens->wz * sens->Bx - sens->wx * sens->Bz;
        double cz =  sens->wx * sens->By - sens->wy * sens->Bx;
        std::cout
          << "t="<<sens->t
          <<" mode="<<(g_mode==CTRL_MODE_WX_B?"WX_B":"B_DOT")
          <<" |B|="<<Bnorm
          <<" ω×b=["<<cx<<","<<cy<<","<<cz<<"]\n";
    }

    // 1) Unpack sensor readings (including timestamp)
    double wx = sens->wx;
    double wy = sens->wy;
    double wz = sens->wz;
    double Bx = sens->Bx;  // in Tesla!
    double By = sens->By;
    double Bz = sens->Bz;
    double t  = sens->t;

    // 2) Compute field norm, bail out if too small
    double Bnorm = std::sqrt(Bx*Bx + By*By + Bz*Bz);
    if (Bnorm < 1e-12) {
        act->mx = act->my = act->mz = 0.0;
        return;
    }

    // 3) Unit‐field vector (Tesla-based)
    double bx = Bx / Bnorm;
    double by = By / Bnorm;
    double bz = Bz / Bnorm;

    // 4) Compute ω×b
    double cx =  wy * bz - wz * by;
    double cy =  wz * bx - wx * bz;
    double cz =  wx * by - wy * bx;

    // 5) Apply control law
    constexpr double k    = 1.0;  // control gain, **unitless here (A·m²/Tesla·rad/s)**
    constexpr double mmax = 0.2;  // dipole‐saturation limit [A·m²]
    double mx = 0.0, my = 0.0, mz = 0.0;

    if (g_mode == CTRL_MODE_B_DOT) {
        // True Bang–bang Ḃ law
        static double Bx_last = 0.0, By_last = 0.0, Bz_last = 0.0, t_last = 0.0;
        double dt = t - t_last;
        if (dt > 0) {
            double Bdotx = (Bx - Bx_last) / dt;
            double Bdoty = (By - By_last) / dt;
            double Bdotz = (Bz - Bz_last) / dt;
            mx = -k * (Bdotx > 0.0 ?  1.0 : (Bdotx < 0.0 ? -1.0 : 0.0));
            my = -k * (Bdoty > 0.0 ?  1.0 : (Bdoty < 0.0 ? -1.0 : 0.0));
            mz = -k * (Bdotz > 0.0 ?  1.0 : (Bdotz < 0.0 ? -1.0 : 0.0));
        }
        // Update memory
        Bx_last = Bx;
        By_last = By;
        Bz_last = Bz;
        t_last  = t;
    } else {
        // Proportional ω×B law
        double factor = k / Bnorm; // ω×b / |B|, still A·m² units
        mx = factor * cx;
        my = factor * cy;
        mz = factor * cz;
    }

    // 6) Saturate each component
    if (mx >  mmax) mx =  mmax; else if (mx < -mmax) mx = -mmax;
    if (my >  mmax) my =  mmax; else if (my < -mmax) my = -mmax;
    if (mz >  mmax) mz =  mmax; else if (mz < -mmax) mz = -mmax;

    // 7) Output commands (A·m²)
    act->mx = mx;
    act->my = my;
    act->mz = mz;
}
