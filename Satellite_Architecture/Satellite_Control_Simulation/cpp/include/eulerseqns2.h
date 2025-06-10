// eulerseqns2.h
#ifndef EULERSEQNS2_H
#define EULERSEQNS2_H

#include <array>

// Global control torque vector L (body frame), to be set by odefun caller
extern std::array<double, 3> L;

// Compute derivative [w_dot; q_dot] for a rigid-body spacecraft
// y: {w1, w2, w3, q1, q2, q3, q4}
// Returns dydt: same structure
std::array<double, 7> eulerseqns2(double t, const std::array<double, 7> &y);

#endif // EULERSEQNS2_H