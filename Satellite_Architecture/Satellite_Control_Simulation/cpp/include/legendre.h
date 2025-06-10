// legendre.hpp

#ifndef LEGENDRE_HPP
#define LEGENDRE_HPP

#include <vector>

// Compute Schmidt semi‐normalized associated Legendre polynomials P_n^m(x)
// and their derivatives dP_n^m/dx at a given x = sin(theta).
//
// Arguments:
//   x      – sin(theta), where theta is colatitude (|x| ≤ 1).
//   P      – pre‐allocated (N+1)×(N+1) matrix; on exit,
//            P[n][m] = P_n^m(x) for 0 ≤ m ≤ n ≤ N.
//   dPdx   – similarly sized; on exit,
//            dPdx[n][m] = ∂P_n^m/∂x.
//
// P and dPdx must be sized so that P.size() = N+1 and P[n].size() = n+1 for each n.
//
// Usage example (inside IGRF code):
//   int N = max_degree;
//   std::vector<std::vector<double>> P(N+1), dPdx(N+1);
//   for (int n = 0; n <= N; ++n) {
//       P[n].resize(n+1);
//       dPdx[n].resize(n+1);
//   }
//   computeLegendre(sin_theta, P, dPdx);
//
// After calling, P[n][m] and dPdx[n][m] are ready for 0 ≤ m ≤ n ≤ N.
//
// Throws if x∉[−1,1] or if P/dPdx have mismatched dimensions.
void computeLegendre(double x,
                     std::vector<std::vector<double>> &P,
                     std::vector<std::vector<double>> &dPdx);

#endif // LEGENDRE_HPP
