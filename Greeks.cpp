#include "Greeks.hpp"

// =============================================================================
//  Greeks.cpp  —  Numerical Greeks computed from the solved PDE grid
// =============================================================================

Greeks GreeksCalculator::compute(const Grid& solved_grid,
                                  const Option& opt, const BSParams& p) {
    Greeks g;
    const int N = solved_grid.N();

    g.S    .resize(N + 1);
    g.delta.resize(N + 1);
    g.gamma.resize(N + 1);
    g.theta.resize(N + 1);
    g.vega .resize(N + 1);
    g.rho  .resize(N + 1);

    for (int i = 0; i <= N; ++i)
        g.S[i] = solved_grid.S(i);

    compute_delta_gamma(g, solved_grid);         // free: spatial FD on V[i]
    compute_theta      (g, solved_grid, p);      // free: PDE identity
    compute_vega       (g, solved_grid, opt, p); // 2 solver calls
    compute_rho        (g, solved_grid, opt, p); // 2 solver calls

    return g;
}

// =============================================================================
//  Delta and Gamma  —  Spatial finite differences
// =============================================================================
//
//  Both Greeks are obtained for free from the already-computed solution V[i].
//  No re-solve needed: we differentiate the solution vector numerically.
//
void GreeksCalculator::compute_delta_gamma(Greeks& g, const Grid& grid) {
    const int    N  = grid.N();
    const double dS = grid.dS();

    // -- Delta: first derivative dV/dS ----------------------------------------
    //  2nd-order centred difference for interior nodes, 1st-order one-sided
    //  at boundaries (where the stencil would reach outside the grid).
    for (int i = 0; i <= N; ++i) {
        if (i == 0)
            g.delta[i] = (grid[1] - grid[0]) / dS;
        else if (i == N)
            g.delta[i] = (grid[N] - grid[N - 1]) / dS;
        else
            g.delta[i] = (grid[i + 1] - grid[i - 1]) / (2.0 * dS);
    }

    // -- Gamma: second derivative d^2V/dS^2 -----------------------------------
    //  2nd-order centred difference for interior nodes.
    //  Boundary nodes are set equal to their nearest interior neighbour
    //  (the payoff is linear at both ends, so Gamma -> 0 there anyway).
    std::vector<double> gamma_interior(N - 1);
    for (int i = 1; i <= N - 1; ++i)
        gamma_interior[i - 1] = (grid[i + 1] - 2.0 * grid[i] + grid[i - 1]) / (dS * dS);

    g.gamma[0] = gamma_interior[0];
    for (int i = 1; i <= N - 1; ++i)
        g.gamma[i] = gamma_interior[i - 1];
    g.gamma[N] = gamma_interior[N - 2];
}

// =============================================================================
//  Theta  —  Via the PDE identity (no re-solve)
// =============================================================================
//
//  The BS PDE can be written as:  dV/dt = -L[V]
//  where the spatial operator is:
//    L[V]_i = 1/2*sigma^2*S_i^2*Gamma_i + r*S_i*Delta_i - r*V_i
//
//  So Theta is read directly from the Greeks we already computed:
//    Theta_i = -[ 1/2*sigma^2*S_i^2*Gamma_i + r*S_i*Delta_i - r*V_i ]
//
//  Theta is returned in years. Divide by 365 to get daily time decay.
//
void GreeksCalculator::compute_theta(Greeks& g, const Grid& grid, const BSParams& p) {
    const int    N      = grid.N();
    const double sigma2 = p.sigma * p.sigma;

    for (int i = 0; i <= N; ++i) {
        const double Si = grid.S(i);
        g.theta[i] = -(0.5 * sigma2 * Si * Si * g.gamma[i]
                     + p.r * Si * g.delta[i]
                     - p.r * grid[i]);
    }
}

// =============================================================================
//  Vega  —  Bump-and-reprice on sigma (central finite difference)
// =============================================================================
//
//  nu_i = [ V(sigma+eps, S_i) - V(sigma-eps, S_i) ] / (2*eps)
//
//  We cannot differentiate sigma from the grid because it is a model parameter,
//  not a grid variable. The only option is to re-solve twice with perturbed
//  sigma and apply a central finite difference.
//
void GreeksCalculator::compute_vega(Greeks& g, const Grid& ref_grid,
                                     const Option& opt, const BSParams& p) {
    const double eps = 1e-3;   // 10 bps in vol
    const int    N   = ref_grid.N();

    Grid up_grid(ref_grid.S_max(), ref_grid.T(), ref_grid.N(), ref_grid.M());
    Grid dn_grid(ref_grid.S_max(), ref_grid.T(), ref_grid.N(), ref_grid.M());

    solver_.solve(up_grid, opt, p.with_sigma(p.sigma + eps));
    solver_.solve(dn_grid, opt, p.with_sigma(p.sigma - eps));

    for (int i = 0; i <= N; ++i)
        g.vega[i] = (up_grid[i] - dn_grid[i]) / (2.0 * eps);
}

// =============================================================================
//  Rho  —  Bump-and-reprice on r (central finite difference)
// =============================================================================
//
//  rho_i = [ V(r+eps, S_i) - V(r-eps, S_i) ] / (2*eps)
//
//  Same structure as compute_vega but bumping the risk-free rate.
//
void GreeksCalculator::compute_rho(Greeks& g, const Grid& ref_grid,
                                    const Option& opt, const BSParams& p) {
    const double eps = 1e-4;   // 1 bp in rates
    const int    N   = ref_grid.N();

    Grid up_grid(ref_grid.S_max(), ref_grid.T(), ref_grid.N(), ref_grid.M());
    Grid dn_grid(ref_grid.S_max(), ref_grid.T(), ref_grid.N(), ref_grid.M());

    solver_.solve(up_grid, opt, p.with_r(p.r + eps));
    solver_.solve(dn_grid, opt, p.with_r(p.r - eps));

    for (int i = 0; i <= N; ++i)
        g.rho[i] = (up_grid[i] - dn_grid[i]) / (2.0 * eps);
}
