#pragma once
#include "Solver.hpp"
#include "Thomas.hpp"

// =============================================================================
//  CrankNicolson.hpp  —  Crank-Nicolson scheme on a uniform S-grid
// =============================================================================
//
//  Solves the Black-Scholes PDE backward in time (t in [0, T]) :
//
//    dV/dt  +  1/2 * sigma^2 * S^2 * d^2V/dS^2  +  r*S * dV/dS  -  r*V  =  0
//
//  Spatial operator on the uniform grid S_i = i*dS :
//
//    L[V]_i = a_i * V_{i-1}  +  b_i * V_i  +  c_i * V_{i+1}
//
//  where a_i, b_i, c_i come from centred finite differences (see .cpp).
//
//  Crank-Nicolson (theta = 1/2): average of implicit and explicit levels
//
//    (V^{n-1}_i - V^n_i) / dt  =  1/2 * (L[V^{n-1}]_i + L[V^n]_i)
//
//  Rearranged into a tridiagonal system  A * V^{n-1} = d  solved by Thomas.
//
//  Key properties:
//    - Second-order accurate in both time O(dt^2) and space O(dS^2)
//    - Unconditionally stable (von Neumann analysis)
//    - O(N) cost per time step via the Thomas algorithm
// =============================================================================

class CrankNicolson : public Solver {
public:
    void        solve(Grid& grid, const Option& opt, const BSParams& p) override;
    const char* name() const override { return "Crank-Nicolson"; }

private:
    // Advance one time step: V^n (stored in grid) -> V^{n-1} (written to grid).
    // t_next is the target time level (closer to 0).
    void step(Grid& grid, const Option& opt, const BSParams& p, double t_next);
};
