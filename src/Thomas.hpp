#pragma once
#include <vector>

// =============================================================================
//  Thomas.hpp  —  Thomas algorithm: O(N) tridiagonal system solver
// =============================================================================
//
//  Solves  A*x = d  where A is a tridiagonal N×N matrix:
//
//    [ b[0]  c[0]                    ] [ x[0]   ]   [ d[0]   ]
//    [ a[1]  b[1]  c[1]              ] [ x[1]   ]   [ d[1]   ]
//    [       a[2]  b[2]  c[2]        ] [ x[2]   ] = [ d[2]   ]
//    [             ...               ] [  ...   ]   [  ...   ]
//    [             a[n-1]  b[n-1]    ] [ x[n-1] ]   [ d[n-1] ]
//
//  Algorithm: Gaussian elimination forward sweep + back substitution.
//  Cost: O(N) — optimal for tridiagonal systems.
//
//  Used by CrankNicolson and ReducedCN at every time step.
// =============================================================================

class ThomasAlgo {
public:
    // Solve the tridiagonal system and return the solution vector x.
    // Precondition: the matrix must be non-singular after elimination
    // (guaranteed for the CN system under standard BS parameters).
    static std::vector<double> solve(
        const std::vector<double>& a,   // sub-diagonal   (a[0] unused)
        const std::vector<double>& b,   // main diagonal
        const std::vector<double>& c,   // super-diagonal (c[n-1] unused)
        const std::vector<double>& d    // right-hand side
    );
};
