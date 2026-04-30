#include "Thomas.hpp"

// =============================================================================
//  Thomas.cpp  —  Thomas algorithm implementation
// =============================================================================
//
//  Two-pass LU factorisation specialised for tridiagonal matrices.
//
//  Pass 1 — Forward sweep: eliminate the sub-diagonal by row operations,
//  producing modified super-diagonal c' and right-hand side d'.
//
//  Pass 2 — Back substitution: solve the resulting upper-bidiagonal system
//  from the last equation upward.
// =============================================================================

std::vector<double> ThomasAlgo::solve(
    const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d
) {
    std::vector<double> c_prime = c;   // modified super-diagonal
    std::vector<double> d_prime = d;   // modified right-hand side

    const int n = static_cast<int>(d_prime.size());
    const int m = n - 1;               // index of last row

    // -- Forward sweep ---------------------------------------------------------
    //  Normalise first row, then eliminate a[i] from each subsequent row.
    c_prime[0] /= b[0];
    d_prime[0] /= b[0];

    for (int i = 1; i < m; ++i) {
        const double denom = b[i] - a[i] * c_prime[i - 1];
        c_prime[i] /= denom;
        d_prime[i]  = (d_prime[i] - a[i] * d_prime[i - 1]) / denom;
    }

    // Last row (no c[m] to normalise)
    d_prime[m] = (d_prime[m] - a[m] * d_prime[m - 1])
               / (b[m] - a[m] * c_prime[m - 1]);

    // -- Back substitution -----------------------------------------------------
    for (int i = m - 1; i >= 0; --i)
        d_prime[i] -= c_prime[i] * d_prime[i + 1];

    return d_prime;   // d_prime now holds the solution x
}
