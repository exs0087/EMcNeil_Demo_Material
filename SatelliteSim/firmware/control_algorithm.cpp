#include <math.h>
#include "control_algorithm.h"

// Global mode and threshold for auto-switch
static ControlMode g_mode = CTRL_MODE_B_DOT;
static const double DESPIN_THRESHOLD = 0.01; // rad/s

extern "C" void control_algorithm_set_mode(ControlMode mode) {
    g_mode = mode;
}

extern "C" void control_algorithm(const SensorReadings* sens,
                                  ActuatorCommands*   act)
{
    // Unpack sensors
    double wx = sens->wx;
    double wy = sens->wy;
    double wz = sens->wz;
    double Bx = sens->Bx;
    double By = sens->By;
    double Bz = sens->Bz;

    // Auto-switch: if in B-dot and spun down, go to ω×B
    double omega_norm = sqrt(wx*wx + wy*wy + wz*wz);
    if (g_mode == CTRL_MODE_B_DOT && omega_norm < DESPIN_THRESHOLD) {
        g_mode = CTRL_MODE_WX_B;
    }

    // Field magnitude
    double Bnorm = sqrt(Bx*Bx + By*By + Bz*Bz);
    if (Bnorm < 1e-12) {
        // No valid field — zero dipoles
        act->mx = act->my = act->mz = 0.0;
        return;
    }

    // Unit field vector
    double bx = Bx / Bnorm;
    double by = By / Bnorm;
    double bz = Bz / Bnorm;

    const double k    = 1.0;  // control gain
    const double mmax = 0.2;  // saturation limit
    double mx=0, my=0, mz=0;

    if (g_mode == CTRL_MODE_B_DOT) {
        // Bang–bang B·dot: m = -k sign(ω×b)
        double cx =  wy*bz - wz*by;
        double cy =  wz*bx - wx*bz;
        double cz =  wx*by - wy*bx;
        mx = -k * (cx > 0 ? 1 : (cx < 0 ? -1 : 0));
        my = -k * (cy > 0 ? 1 : (cy < 0 ? -1 : 0));
        mz = -k * (cz > 0 ? 1 : (cz < 0 ? -1 : 0));
    }
    else {
        // ω×B law (despin): m = (k/Bnorm)*(ω×b)
        double cx =  wy*bz - wz*by;
        double cy =  wz*bx - wx*bz;
        double cz =  wx*by - wy*bx;
        double factor = k / Bnorm;
        mx = factor * cx;
        my = factor * cy;
        mz = factor * cz;
    }

    // Saturation
    if (mx >  mmax) mx =  mmax; else if (mx < -mmax) mx = -mmax;
    if (my >  mmax) my =  mmax; else if (my < -mmax) my = -mmax;
    if (mz >  mmax) mz =  mmax; else if (mz < -mmax) mz = -mmax;

    act->mx = mx;
    act->my = my;
    act->mz = mz;
}
