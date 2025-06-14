#include "eulerseqns2.h"
#include "magnetic_field.h"
#include "Vec3.h"
#include "Vec4.h"
#include <limits>

// ----------------------------------------------------------------------------
// Globals to collect one actuator command per unique top‐level time‐step
static std::vector<ActuatorCommands> g_act_hist;
static double                       g_last_recorded_t = std::numeric_limits<double>::quiet_NaN();

// ----------------------------------------------------------------------------
Vec7 eulerseqns2(double t, const Vec7& y) {
    // 1) Unpack state
    Vec3 w{ y[0], y[1], y[2] };           // body rates [rad/s]
    Vec4 q{ y[3], y[4], y[5], y[6] };      // quaternion (scalar last)
    q.normalize();

    // 2) Get body‐frame magnetic field from the sim HAL
    Vec3 B_body = getBodyFieldAt(t, q);

    // 3) Build sensor reading struct
    SensorReadings sens;
    sens.wx = w.x;  sens.wy = w.y;  sens.wz = w.z;
    sens.Bx = B_body.x;  sens.By = B_body.y;  sens.Bz = B_body.z;

    // 4) Call your embedded control algorithm
    ActuatorCommands act{0,0,0};
    HAL_compute_control(&sens, &act);

    // 5) Record exactly one command per unique t (first call at that timestamp)
    if (t != g_last_recorded_t) {
        g_act_hist.push_back(act);
        g_last_recorded_t = t;
    }

    // 6) Compute control torque: τ = m × B
    Vec3 m{ act.mx, act.my, act.mz };
    Vec3 torque = m.cross(B_body);

    // 7) Rigid‐body Euler equations + quaternion kinematics
    const double mass  = 10.0, hx = 0.34, wy = 0.20, dz = 0.10;
    const double Jx    = mass/12*(wy*wy + dz*dz);
    const double Jy    = mass/12*(hx*hx + dz*dz);
    const double Jz    = mass/12*(hx*hx + wy*wy);
    const double invJx = 1.0/Jx, invJy = 1.0/Jy, invJz = 1.0/Jz;
    const double c1    = (Jy - Jz)*invJx;
    const double c2    = (Jz - Jx)*invJy;
    const double c3    = (Jx - Jy)*invJz;

    Vec7 dydt;
    dydt[0] = c1*w.y*w.z + torque.x*invJx;
    dydt[1] = c2*w.z*w.x + torque.y*invJy;
    dydt[2] = c3*w.x*w.y + torque.z*invJz;

    // Quaternion derivative: q̇ = ½ Ω(ω) q
    double ωx = w.x, ωy = w.y, ωz = w.z;
    dydt[3] =  0.5*(  q.w*ωx + q.y*ωz - q.z*ωy );
    dydt[4] =  0.5*(  q.w*ωy + q.z*ωx - q.x*ωz );
    dydt[5] =  0.5*(  q.w*ωz + q.x*ωy - q.y*ωx );
    dydt[6] =  0.5*( -q.x*ωx - q.y*ωy - q.z*ωz );

    return dydt;
}

// ----------------------------------------------------------------------------
const std::vector<ActuatorCommands>& getActHistory() {
    return g_act_hist;
}

// ----------------------------------------------------------------------------
void resetActHistory() {
    g_act_hist.clear();
    g_last_recorded_t = std::numeric_limits<double>::quiet_NaN();
}
