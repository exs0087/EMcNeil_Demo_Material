#pragma once
#include <array>

/// 7-element state: [ω₁, ω₂, ω₃, q₁, q₂, q₃, q₄]
using Vec7 = std::array<double,7>;
/// 3-element control torque vector: [L₁, L₂, L₃]
using Vec3 = std::array<double,3>;

/**
 * Compute the rigid-body Euler rotation plus quaternion kinematics.
 * This is a direct port of eulerseqns2.m:
 *   ω̇ = I⁻¹ (L - ω × (I ω))
 *   q̇ = ½ Ω(ω) q
 *
 * @param t  Current time (unused)
 * @param y  State vector [ω₁,ω₂,ω₃,q₁,q₂,q₃,q₄]
 * @param L  Control torque vector [L₁,L₂,L₃] (N·m)
 * @returns  dydt = [ω̇₁,ω̇₂,ω̇₃,q̇₁,q̇₂,q̇₃,q̇₄]
 */
Vec7 eulerDynamics(double t, const Vec7 &y, const Vec3 &L);