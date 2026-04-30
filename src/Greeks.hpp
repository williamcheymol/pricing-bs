#pragma once
#include <vector>
#include "Grid.hpp"
#include "Solver.hpp"
#include "Option.hpp"

// =============================================================================
//  Greeks.hpp  —  Finite-difference Greeks from the solved PDE grid
// =============================================================================
//
//  Greeks measure the sensitivity of the option price V to model inputs.
//  Three different methods are used depending on what the Greek differentiates:
//
//  +---------+-----------------------------------+-------------------------------+
//  | Greek   | Method                            | Cost                          |
//  +---------+-----------------------------------+-------------------------------+
//  | Delta   | Central FD in S on solved grid    | Free (reuse V[i])             |
//  | Gamma   | Central FD in S on solved grid    | Free (reuse V[i])             |
//  | Theta   | PDE identity: Theta = -L[V]       | Free (reuse Delta, Gamma)     |
//  | Vega    | Bump-and-reprice on sigma (+-eps) | 2 full solver calls           |
//  | Rho     | Bump-and-reprice on r     (+-eps) | 2 full solver calls           |
//  +---------+-----------------------------------+-------------------------------+
//
//  Total cost per option: 5 solver calls (1 main + 2 Vega + 2 Rho).
//
//  Validation: compare each FD Greek against BSAnalytical (BSFormula.hpp).
// =============================================================================

struct Greeks {
    std::vector<double> S;      // spot grid (reference axis)
    std::vector<double> delta;  // Delta = dV/dS
    std::vector<double> gamma;  // Gamma = d^2V/dS^2
    std::vector<double> theta;  // Theta = dV/dt  (per year; divide by 365 for daily)
    std::vector<double> vega;   // Vega  = dV/d(sigma)
    std::vector<double> rho;    // Rho   = dV/dr
};

class GreeksCalculator {
public:
    // The solver reference is kept for bump-and-reprice (Vega, Rho).
    // Any Solver implementation (CrankNicolson, ReducedCN, ...) works here.
    explicit GreeksCalculator(Solver& solver) : solver_(solver) {}

    // Compute all Greeks on the full spatial grid and return them.
    // solved_grid must already contain V(t=0, S_i) from a prior solve() call.
    Greeks compute(const Grid& solved_grid, const Option& opt, const BSParams& p);

private:
    Solver& solver_;

    // -- Delta & Gamma : spatial finite differences ----------------------------
    //
    //  Interior nodes (2nd-order centred):
    //    Delta_i = (V_{i+1} - V_{i-1}) / (2*dS)
    //    Gamma_i = (V_{i+1} - 2*V_i + V_{i-1}) / dS^2
    //
    //  Boundary nodes (1st-order one-sided for Delta, extension for Gamma):
    //    Delta_0 = (V_1 - V_0) / dS          Gamma_0 = Gamma_1
    //    Delta_N = (V_N - V_{N-1}) / dS      Gamma_N = Gamma_{N-1}
    //
    void compute_delta_gamma(Greeks& g, const Grid& grid);

    // -- Theta : PDE identity (no extra solve needed) --------------------------
    //
    //  From the BS PDE:  dV/dt = -L[V]
    //  with  L[V]_i = 1/2*sigma^2*S_i^2*Gamma_i + r*S_i*Delta_i - r*V_i
    //
    //  Therefore:  Theta_i = -[ 1/2*sigma^2*S_i^2*Gamma_i
    //                           + r*S_i*Delta_i - r*V_i ]
    //
    //  More accurate than finite-differencing in time: no extra storage needed
    //  and exploits the exact PDE relationship between the spatial Greeks.
    //
    void compute_theta(Greeks& g, const Grid& grid, const BSParams& p);

    // -- Vega : bump-and-reprice on sigma (central FD) -------------------------
    //
    //  nu_i = [ V(sigma+eps, S_i) - V(sigma-eps, S_i) ] / (2*eps)
    //
    //  eps = 10 bps in vol: large enough to avoid floating-point cancellation,
    //  small enough for the finite-difference approximation to be accurate.
    //
    void compute_vega(Greeks& g, const Grid& ref_grid,
                      const Option& opt, const BSParams& p);

    // -- Rho : bump-and-reprice on r (central FD) ------------------------------
    //
    //  rho_i = [ V(r+eps, S_i) - V(r-eps, S_i) ] / (2*eps)
    //
    //  eps = 1 bp in rates.
    //
    void compute_rho(Greeks& g, const Grid& ref_grid,
                     const Option& opt, const BSParams& p);
};
