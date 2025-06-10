// eulerseqns2.cpp
#include "eulerseqns2.h"
#include <cmath>

// Definition of global L; must be set in the simulation before each derivative call
std::array<double, 3> L = {0.0, 0.0, 0.0};

std::array<double, 7> eulerseqns2(double t, const std::array<double, 7> &y) {
    // Mass properties (rectangular prism)
    double mass = 10.0;      // kg
    double height = 0.34;    // m (x dimension)
    double width  = 0.20;    // m (y dimension)
    double depth  = 0.10;    // m (z dimension)

    // Principal moments of inertia
    double Jx = mass / 12.0 * (width*width + depth*depth);
    double Jy = mass / 12.0 * (height*height + depth*depth);
    double Jz = mass / 12.0 * (height*height + width*width);

    // Inertia coefficients for Euler equations
    double c1 = (Jy - Jz) / Jx;
    double c2 = (Jz - Jx) / Jy;
    double c3 = (Jx - Jy) / Jz;

    // Unpack state
    double w1 = y[0];
    double w2 = y[1];
    double w3 = y[2];
    std::array<double,3> w = { w1, w2, w3 };
    double q1 = y[3];
    double q2 = y[4];
    double q3 = y[5];
    double q4 = y[6];

    // Construct xi matrix (4x3) to compute quaternion derivative
    // xi * w = q_dot * 2
    // xi = [ [ q4, -q3,  q2 ];
    //        [ q3,  q4, -q1 ];
    //        [-q2,  q1,  q4 ];
    //        [-q1, -q2, -q3 ] ]
    double xi[4][3] = {
        {  q4, -q3,  q2 },
        {  q3,  q4, -q1 },
        { -q2,  q1,  q4 },
        { -q1, -q2, -q3 }
    };

    // Compute q_dot = 0.5 * xi * w
    std::array<double,4> qd;
    for (int i = 0; i < 4; ++i) {
        qd[i] = 0.5 * ( xi[i][0] * w1 + xi[i][1] * w2 + xi[i][2] * w3 );
    }

    // Compute angular acceleration from Euler's equations: w_dot = [c1 w2 w3; c2 w3 w1; c3 w1 w2] + [Lx/Jx; Ly/Jy; Lz/Jz]
    std::array<double,3> wdot;
    wdot[0] = c1 * w2 * w3 + L[0] / Jx;
    wdot[1] = c2 * w3 * w1 + L[1] / Jy;
    wdot[2] = c3 * w1 * w2 + L[2] / Jz;

    // Assemble output derivative vector
    std::array<double,7> dydt;
    dydt[0] = wdot[0];
    dydt[1] = wdot[1];
    dydt[2] = wdot[2];
    dydt[3] = qd[0];
    dydt[4] = qd[1];
    dydt[5] = qd[2];
    dydt[6] = qd[3];

    return dydt;
}