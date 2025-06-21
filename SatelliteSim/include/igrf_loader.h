#pragma once

#include <vector>

/// Load and interpolate the IGRF g/h coefficient vector for a given decimal year.
/// Reads epoch years from data/years.csv and coefficient vectors from data/gh.csv.
/// @param decimalYear  Year in decimal (e.g. 2025.5)
/// @returns            Interpolated [g₁⁰, g₁¹, h₁¹, …, gₙᵐ, hₙᵐ] vector
std::vector<double> loadIGRFCoeffs(double decimalYear);
