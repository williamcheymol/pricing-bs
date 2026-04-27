#include "ReducedCN.hpp"
#include "Thomas.hpp"
#include <cmath>
#include <algorithm>

// Portable clamp for integer indices (std::clamp requires C++17 which some
// toolchains don't enable by default).
static inline int iclamp(int v, int lo, int hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

// =============================================================================
//  ReducedCN::solve  —  Full backward solve on the log-uniform grid
// =============================================================================

void ReducedCN::solve(Grid& grid, const Option& opt, const BSParams& p) {
    const int    N     = grid.N();
    const double S_max = grid.S_max();
    const double T     = p.T;
    const int    M     = grid.M();
    const double dt    = T / M;
    const double dS    = S_max / N;

    // -- Build the log-uniform grid  x in [x_min, x_max]  with N steps --------
    //
    //  x_min = ln(dS) = ln(S_max / N)   effective lower boundary  (S = dS)
    //  x_max = ln(S_max)                upper boundary             (S = S_max)
    //
    //  Why x_min = ln(dS) and not ln(0)?
    //  S = 0 maps to x = -inf and cannot be represented on a finite grid.
    //  We use S_1 = dS as the effective lower node.  For N = 1000, dS = 0.30,
    //  so the approximation bc_lower(dS) ≈ bc_lower(0) is excellent.
    //
    //  The S-values on this log-grid are geometrically spaced:
    //    S_i^{log} = dS * N^{i/N}   (ratio dS/S stays constant across nodes)
    const double x_min = std::log(dS);
    const double x_max = std::log(S_max);
    const double dx    = (x_max - x_min) / N;

    // -- Terminal condition: V(T, x_i) = payoff(e^{x_i}) ----------------------
    std::vector<double> V(N + 1);
    for (int i = 0; i <= N; ++i)
        V[i] = opt.payoff(std::exp(x_min + i * dx));

    // -- Backward time march from t = T to t = 0 ------------------------------
    for (int m = M; m >= 1; --m) {
        double t_next = (m - 1) * dt;
        step(V, opt, p, x_min, dx, N, dt, S_max, t_next);
    }

    // -- Project solution onto the uniform S-grid ------------------------------
    //
    //  The output grid uses S_j = j*dS (j = 0, ..., N).
    //  For each S_j, compute x_j = ln(S_j) and interpolate linearly within
    //  the log-grid {x_i = x_min + i*dx}.
    //
    //  Special cases:
    //    j = 0 : S = 0 -> exact analytical BC (ln(0) is undefined)
    //    j = 1 : x = ln(dS) = x_min -> maps exactly to V[0]
    //    j = N : x = ln(S_max) = x_max -> maps exactly to V[N]
    grid[0] = opt.bc_lower(0.0, T, p.r);
    for (int j = 1; j <= N; ++j) {
        double xj    = std::log(static_cast<double>(j) * dS);
        double fi    = (xj - x_min) / dx;   // fractional position in log-grid
        int    il    = iclamp(static_cast<int>(fi), 0, N - 1);
        double alpha = fi - il;             // interpolation weight in [0, 1)
        grid[j] = (1.0 - alpha) * V[il] + alpha * V[il + 1];
    }
}

// =============================================================================
//  ReducedCN::step  —  One time step on the log-grid
// =============================================================================
//
//  Reduced spatial operator (CONSTANT coefficients):
//    L[V]_i = a * V_{i-1}  +  b * V_i  +  c * V_{i+1}
//
//  where kappa = r - 1/2*sigma^2  and:
//    a =  dt/4 * (sigma^2/dx^2  -  kappa/dx)
//    b = -dt/2 * (sigma^2/dx^2  +  r       )
//    c =  dt/4 * (sigma^2/dx^2  +  kappa/dx)
//
//  Same CN structure as CrankNicolson::step(), but a, b, c are constant
//  across all nodes — no S_i dependence, no coefficient growth for large S.
//
//  Tridiagonal system (unknowns at level n-1, knowns at level n):
//    -a * V^{n-1}_{i-1} + (1-b) * V^{n-1}_i - c * V^{n-1}_{i+1}
//    =  a * V^n_{i-1}   + (1+b) * V^n_i      + c * V^n_{i+1}
//
//  Boundary correction (same logic as CrankNicolson::step):
//    k=0   : rhs[0]   -= lower[0]   * bc_lo
//    k=n-1 : rhs[n-1] -= upper[n-1] * bc_hi
//
// =============================================================================
void ReducedCN::step(std::vector<double>& V,
                     const Option&        opt,
                     const BSParams&      p,
                     double               /*x_min*/,
                     double               dx,
                     int                  Nx,
                     double               dt,
                     double               S_max,
                     double               t_next) {
    const int    n      = Nx - 1;
    const double sigma2 = p.sigma * p.sigma;
    const double kappa  = p.r - 0.5 * sigma2;   // log-space drift

    // -- Constant operator coefficients ----------------------------------------
    const double a = 0.25 * dt * (sigma2 / (dx * dx) - kappa / dx);
    const double b = -0.5 * dt * (sigma2 / (dx * dx) + p.r);
    const double c = 0.25 * dt * (sigma2 / (dx * dx) + kappa / dx);

    // -- Assemble tridiagonal system --------------------------------------------
    //  lower and upper are uniform (constant a, c) — unlike the S-space CN.
    std::vector<double> lower(n, -a);
    std::vector<double> diag (n, 1.0 - b);
    std::vector<double> upper(n, -c);
    std::vector<double> rhs  (n);

    for (int k = 0; k < n; ++k) {
        int i  = k + 1;
        rhs[k] = a * V[i - 1] + (1.0 + b) * V[i] + c * V[i + 1];
    }

    // -- Apply boundary conditions at t_next -----------------------------------
    const double T    = p.T;
    double bc_lo = opt.bc_lower(t_next, T, p.r);
    double bc_hi = opt.bc_upper(S_max, t_next, T, p.r);

    V[0]  = bc_lo;
    V[Nx] = bc_hi;

    rhs[0]     -= lower[0]     * bc_lo;
    rhs[n - 1] -= upper[n - 1] * bc_hi;

    // -- Solve in O(N) via Thomas algorithm ------------------------------------
    std::vector<double> sol = ThomasAlgo::solve(lower, diag, upper, rhs);
    for (int k = 0; k < n; ++k)
        V[k + 1] = sol[k];
}
