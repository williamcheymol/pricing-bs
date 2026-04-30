#pragma once
#include "Grid.hpp"
#include "Option.hpp"

// =============================================================================
//  Solver.hpp  —  Abstract solver interface + Black-Scholes parameters
// =============================================================================

// Black-Scholes model parameters.
// Passed by value so solvers can't accidentally mutate the caller's state.
// The helper methods return modified copies for bump-and-reprice Greeks.
struct BSParams {
    double r;      // risk-free rate (continuous, annualised)
    double sigma;  // volatility (annualised)
    double T;      // maturity (in years)

    // Return a copy with sigma bumped — used by GreeksCalculator for Vega.
    BSParams with_sigma(double s) const { return {r, s, T}; }

    // Return a copy with r bumped — used by GreeksCalculator for Rho.
    BSParams with_r(double x) const { return {x, sigma, T}; }
};

// Abstract interface for all numerical PDE schemes.
//
// Contract: after solve(), grid[i] = V(t=0, S_i) for all i.
//
// The interface decouples the Greeks calculator and main driver from any
// specific scheme — CrankNicolson and ReducedCN are interchangeable.
class Solver {
public:
    virtual ~Solver() = default;

    virtual void        solve(Grid& grid, const Option& opt, const BSParams& p) = 0;
    virtual const char* name() const = 0;
};
