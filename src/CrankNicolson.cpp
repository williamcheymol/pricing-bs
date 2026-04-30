#include "CrankNicolson.hpp"
#include <vector>

// =============================================================================
//  CrankNicolson.cpp  —  Crank-Nicolson scheme implementation
// =============================================================================

// Backward time march from t = T down to t = 0.
// Each call to step() solves one tridiagonal system and overwrites the grid.
void CrankNicolson::solve(Grid& grid, const Option& opt, const BSParams& p) {
    grid.initialize(opt);   // set V(T, S_i) = payoff(S_i)

    for (int m = grid.M(); m >= 1; --m) {
        double t_next = (m - 1) * grid.dt();   // target time level (smaller)
        step(grid, opt, p, t_next);
    }
    // On exit: grid[i] = V(t=0, S_i)
}

// Advance one time step using the Crank-Nicolson scheme.
void CrankNicolson::step(Grid& grid, const Option& opt, const BSParams& p, double t_next) {
    const int    N  = grid.N();
    const double dS = grid.dS();
    const double dt = grid.dt();
    const double r  = p.r;

    // Interior nodes i = 1, ..., N-1  ->  n = N-1 unknowns
    const int n = N - 1;
    std::vector<double> lower(n), diag(n), upper(n), rhs(n);

    for (int k = 0; k < n; ++k) {
        const int    i  = k + 1;        // true grid index
        const double Si = grid.S(i);    // S_min + i * dS  (correct for down barriers)

        // -- Spatial operator coefficients L[V]_i = a_i*V_{i-1} + b_i*V_i + c_i*V_{i+1}
        //
        //  Centred FD discretisation of the BS operator:
        //    diffusion  : 1/2*sigma^2*S^2 * (V_{i+1} - 2V_i + V_{i-1}) / dS^2
        //    convection :           r*S   * (V_{i+1} - V_{i-1})         / (2*dS)
        //    discounting:                  - r * V_i
        //
        //  Each coefficient already absorbs the dt/2 factor from the CN average.

        double a_i = 0.25 * dt * ((p.sigma*p.sigma*Si*Si)/(dS*dS) - (r*Si)/dS);
        double b_i = -0.5 * dt * ((p.sigma*p.sigma*Si*Si)/(dS*dS) + r);
        double c_i = 0.25 * dt * ((p.sigma*p.sigma*Si*Si)/(dS*dS) + (r*Si)/dS);

        // -- Tridiagonal system assembly
        //
        //  CN scheme:  (V^{n-1}_i - V^n_i) / dt = 1/2 * (L[V^{n-1}]_i + L[V^n]_i)
        //
        //  Unknowns (level n-1) on the left, knowns (level n) on the right:
        //    -a_i * V^{n-1}_{i-1} + (1-b_i) * V^{n-1}_i - c_i * V^{n-1}_{i+1}
        //    = a_i * V^n_{i-1}   + (1+b_i) * V^n_i      + c_i * V^n_{i+1}

        lower[k] = -a_i;
        diag [k] = 1.0 - b_i;
        upper[k] = -c_i;
        rhs  [k] = a_i * grid[i-1] + (1.0 + b_i) * grid[i] + c_i * grid[i+1];
    }

    // -- Update boundary conditions to the new time level
    grid.apply_bc(opt, t_next, r);

    // -- Incorporate BCs into the RHS
    //
    //  grid[0] and grid[N] are now fixed at t_next.
    //  The equations for i=1 (k=0) and i=N-1 (k=n-1) reference these
    //  boundary values, so we move them to the right-hand side.

    rhs[0]     -= lower[0]     * grid[0];   // lower boundary contribution
    rhs[n - 1] -= upper[n - 1] * grid[N];   // upper boundary contribution

    // -- Solve the tridiagonal system in O(N)
    std::vector<double> sol = ThomasAlgo::solve(lower, diag, upper, rhs);

    // Write interior solution back to the grid
    for (int k = 0; k < n; ++k)
        grid[k + 1] = sol[k];

    // Early-exercise projection (American options) --------------------
    //
    //  After the Thomas solve, grid[i] holds the "continuation value" V*_i —
    //  what the option is worth if the holder does NOT exercise today.
    //
    //  For American options, the holder compares this with the intrinsic value
    //  (the payoff from exercising immediately). The option is worth the max:
    //
    //      V_i = max(V*_i,  opt.intrinsic(S_i))
    //
    //  For European options, opt.intrinsic() always returns 0.0 so this step
    //  is a no-op — implement it unconditionally and it works for everyone.
    //
    for (int i = 0; i <= N; ++i)
        grid[i] = std::max(grid[i], opt.intrinsic(grid.S(i)));
}
