// src/eulerseqns2.cpp

#include "eulerseqns2.h"
#include "iss_orbit.h"
#include "magnetic_field.h"
#include <cmath>
#include <vector>

// Local cross‐product helper, same as in Vec3.h
static Vec3 cross(const Vec3 &a, const Vec3 &b) {
    return Vec3{
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

// History buffers (populated by HAL_compute_control call)
static std::vector<ActuatorCommands> g_dipole_hist;
static std::vector<Vec3>             g_torque_hist;

void clearActHistory() {
    g_dipole_hist.clear();
    g_torque_hist.clear();
}

const std::vector<ActuatorCommands>& getActHistory() {
    return g_dipole_hist;
}

const std::vector<Vec3>& getTorqueHistory() {
    return g_torque_hist;
}

// Earth & spacecraft inertia
static constexpr double mass   = 10.0;    // kg
static constexpr double height = 0.34;    // m (x)
static constexpr double width  = 0.20;    // m (y)
static constexpr double depth  = 0.10;    // m (z)

// Inertia tensor (rectangular prism, kg*m^2)
static constexpr double Jx = mass/12.0 * (width*width + depth*depth);
static constexpr double Jy = mass/12.0 * (height*height + depth*depth);
static constexpr double Jz = mass/12.0 * (height*height + width*width);


Vec7 eulerseqns2(double t, const Vec7& y) {
    Vec7 dydt;

    // --- 1) Unpack attitude state ---
    double wx = y[0], wy = y[1], wz = y[2];
    double q1 = y[3], q2 = y[4], q3 = y[5], q4 = y[6];  // scalar-last

    // --- 2) Magnetic field in body frame ---
    // 2a) ECI pos
    Vec3 pos_eci = issOrbitECI(t);

    // 2b) Geodetic coords are handled inside getBodyFieldAt, which returns B-field in ECI frame (Tesla)
    double lat, lon, alt_km;
    issOrbitGeod(t, lat, lon, alt_km);
    Vec3 Beci_T = getBodyFieldAt(t, Vec4{lat, lon, alt_km, 0.0});

    // 2c) Quaternion → body rotation (still in Tesla)
    // build DCM from q = [qx,qy,qz,qw]
    double qw=q4, qx=q1, qy=q2, qz=q3;
    double r11 = 1-2*(qy*qy+qz*qz), r12 = 2*(qx*qy - qz*qw), r13 = 2*(qx*qz + qy*qw);
    double r21 = 2*(qx*qy + qz*qw), r22 = 1-2*(qx*qx+qz*qz), r23 = 2*(qy*qz - qx*qw);
    double r31 = 2*(qx*qz - qy*qw), r32 = 2*(qy*qz + qx*qw), r33 = 1-2*(qx*qx+qy*qy);
    Vec3 Bbody_T{
        r11*Beci_T.x + r12*Beci_T.y + r13*Beci_T.z,
        r21*Beci_T.x + r22*Beci_T.y + r23*Beci_T.z,
        r31*Beci_T.x + r32*Beci_T.y + r33*Beci_T.z
    };

    // --- 3) Call into HAL to get dipole commands ---
    // Bfield in body frame is in Tesla!
    SensorReadings sens{ wx, wy, wz, Bbody_T.x, Bbody_T.y, Bbody_T.z, t };
    ActuatorCommands act;
    HAL_compute_control(&sens, &act);

    // record dipole
    g_dipole_hist.push_back(act);

    // --- 4) Compute torque L = m × B ---
    // Actuator dipoles assumed SI-units [A·m^2]; B in [T] → torque [N·m]
    Vec3 mvec{ act.mx, act.my, act.mz };
    Vec3 L = cross(mvec, Bbody_T);
    g_torque_hist.push_back(L);

    // --- 5) Attitude dynamics ---
    // Euler’s eqns: I·ω̇ = L – ω×(I·ω)
    dydt[0] = ( L.x - (Jy-Jz)*wy*wz ) / Jx;
    dydt[1] = ( L.y - (Jz-Jx)*wx*wz ) / Jy;
    dydt[2] = ( L.z - (Jx-Jy)*wx*wy ) / Jz;

    // quaternion kinematics (scalar-last)
    dydt[3] =  0.5*(  qw*wx + qz*wy - qy*wz );
    dydt[4] =  0.5*(  qw*wy - qz*wx + qx*wz );
    dydt[5] =  0.5*(  qw*wz + qy*wx - qx*wy );
    dydt[6] =  0.5*( -qx*wx - qy*wy - qz*wz );

    return dydt;
}
