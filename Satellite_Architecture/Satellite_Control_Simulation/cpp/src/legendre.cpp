// legendre.cpp

#include "legendre.hpp"
#include <cmath>
#include <stdexcept>

// Compute Schmidt semi‐normalized associated Legendre functions P_n^m(x) and their
// derivatives ∂P_n^m/∂x at a given x = sin(theta).
//
// Inputs:
//   x        – sin(theta) (double in [−1, 1])
//   P        – a (N+1)×(N+1) matrix (vector of vectors) to be filled with:
//                P[n][m] = P_n^m(x),  for 0 ≤ m ≤ n ≤ N.
//   dPdx     – a similarly‐sized matrix to be filled with ∂P_n^m/∂x.
//
// Both P and dPdx must be pre‐allocated to size (N+1)×(N+1); only entries with
// m ≤ n are written.
//
// The Schmidt semi‐normalization is defined such that
//   ∫_{0}^{π} [P_n^m(cosθ)]^2 sinθ dθ = 1.
// In terms of x = sinθ, this matches the "sch" flag in MATLAB's legendre(..., "sch").
//
// Uses standard recurrence:
//
//   P_0^0(x) = 1
//   P_1^0(x) =  √(1−x²)
//   P_n^0(x) = ( (2n−1) x P_{n−1}^0(x) − (n−1) P_{n−2}^0(x) ) / n
//
//   P_n^n(x) = √(1−x²) * √( (2n−1)/(2n) ) * P_{n−1}^{n−1}(x)
//   P_n^{n−1}(x) = x * √(2n−1) * P_{n−1}^{n−1}(x)
//
//   For 0 < m < n:
//   P_n^m(x) = √( (2n−1)/(n−m)·(n+m) ) * ( x * P_{n−1}^m(x) − √( ((n−1)² − m²)/(2n−3)/(n+m−1)/(n−m−1) ) * P_{n−2}^m(x) )
//
// Derivatives use:
//   ∂P_n^m/∂x = (n x P_n^m(x) − (n+m) P_{n−1}^m(x)) / (x² − 1)
//
// Reference: "Numerical Recipes" and standard geodesy texts for Schmidt semi‐normalized recurrence.

void computeLegendre(double x,
                     std::vector<std::vector<double>> &P,
                     std::vector<std::vector<double>> &dPdx)
{
    // Ensure x in valid range
    if (x < -1.0 || x > 1.0) {
        throw std::domain_error("computeLegendre: x = sin(theta) must be in [-1,1]");
    }

    size_t N = P.size() - 1;
    if (P.size() != dPdx.size()) {
        throw std::invalid_argument("computeLegendre: P and dPdx must have same dimensions");
    }
    for (size_t n = 0; n <= N; ++n) {
        if (P[n].size() != n+1 || dPdx[n].size() != n+1) {
            throw std::invalid_argument("computeLegendre: P and dPdx rows must be length n+1");
        }
    }

    // Precompute sqrt(1 - x^2)
    double one_minus_x2 = 1.0 - x*x;
    double sx = (one_minus_x2 <= 0.0 ? 0.0 : std::sqrt(one_minus_x2));

    // Base case: P_0^0 = 1
    P[0][0] = 1.0;
    dPdx[0][0] = 0.0;

    if (N == 0) return;

    // n = 1, m = 0 or 1
    // P_1^0(x) = sqrt(1 - x^2)
    P[1][0] = sx;
    dPdx[1][0] = -x / sx; // derivative of sqrt(1-x^2) = -x / sqrt(1-x^2)

    // P_1^1(x) = x * sqrt(3) * P_0^0 = sqrt(3)*(x)
    double f = std::sqrt(3.0);
    P[1][1] = f * x;
    dPdx[1][1] = f;

    // Now iterate n = 2..N
    for (size_t n = 2; n <= N; ++n) {
        // First compute m = 0:
        // P_n^0(x) = ( (2n-1) x P_{n-1}^0(x) - (n-1) P_{n-2}^0(x) ) / n
        double a0 = (2.0*n - 1.0) / double(n);
        double b0 = double(n - 1) / double(n);
        P[n][0] = a0 * x * P[n-1][0] - b0 * P[n-2][0];
        // Derivative via: dP_n^0/dx = (n x P_n^0 - n P_{n-1}^0) / (x^2 - 1)
        if (one_minus_x2 > 0.0) {
            dPdx[n][0] = ( double(n)*x*P[n][0] - double(n)*P[n-1][0] ) / (x*x - 1.0);
        } else {
            // At poles, use recursion on derivative of P_n^0:
            dPdx[n][0] = 0.0;
        }

        // Then m = n:
        // P_n^n(x) = sx * sqrt((2n-1)/(2n)) * P_{n-1}^{n-1}(x)
        double cnn = std::sqrt((2.0*n - 1.0) / (2.0*n));
        P[n][n] = sx * cnn * P[n-1][n-1];
        // Derivative: dP_n^n/dx = cnn * [ sqrt(1-x^2)' * P_{n-1}^{n-1} + sx * dP_{n-1}^{n-1}/dx ]
        double dsx_dx = (one_minus_x2 <= 0.0 ? 0.0 : -x / sx);
        if (n == 1) {
            dPdx[n][n] = cnn * dsx_dx * P[n-1][n-1];
        } else {
            dPdx[n][n] = cnn * ( dsx_dx * P[n-1][n-1] + sx * dPdx[n-1][n-1] );
        }

        // And m = n-1:
        // P_n^{n-1}(x) = x * sqrt(2n-1) * P_{n-1}^{n-1}(x)
        double cn1 = std::sqrt(2.0*n - 1.0);
        P[n][n-1] = cn1 * x * P[n-1][n-1];
        // Derivative: dP_n^{n-1}/dx = cn1 * [ P_{n-1}^{n-1} + x * dP_{n-1}^{n-1}/dx ]
        dPdx[n][n-1] = cn1 * ( P[n-1][n-1] + x * dPdx[n-1][n-1] );

        // Finally, for 1 ≤ m ≤ n-2:
        for (size_t m = 1; m <= n - 2; ++m) {
            // Coefficients:
            // α = sqrt( (2n-1)/(n-m)/(n+m) )
            // β = sqrt( ((n-1)^2 - m^2)/(2n-3)/(n+m-1)/(n-m-1) )
            double num = 2.0*n - 1.0;
            double den = double(n - m)*double(n + m);
            double alpha = std::sqrt(num / den);

            double num2 = double((n - 1)*(n - 1)) - double(m*m);
            double den2 = (2.0*n - 3.0) * double(n + m - 1) * double(n - m - 1);
            double beta = (den2 > 0.0 ? std::sqrt(num2 / den2) : 0.0);

            // Recurrence:
            // P_n^m(x) = α [ x * P_{n-1}^m(x) - β * P_{n-2}^m(x) ]
            double term1 = x * P[n-1][m];
            double term2 = beta * P[n-2][m];
            P[n][m] = alpha * (term1 - term2);

            // Derivative: use exact formula:
            // dP_n^m/dx = (n x P_n^m(x) - (n+m) P_{n-1}^m(x)) / (x^2 - 1)
            if (one_minus_x2 > 0.0) {
                dPdx[n][m] = ( double(n)*x*P[n][m] - double(n + m)*P[n-1][m] ) / (x*x - 1.0);
            } else {
                // At poles (|x|=1), derivative is zero for m>0
                dPdx[n][m] = 0.0;
            }
        }
    }
}
