#pragma once
#include <vector>
#include "Option.hpp"

// =============================================================================
//  Grid.hpp  —  Uniform spatial / temporal discretisation grid
// =============================================================================
//
//  Uniform grid over (S, t) :
//    S_i = S_min + i * dS,   i = 0, 1, ..., N     (N+1 spatial nodes)
//    t_m = m * dt,           m = 0, 1, ..., M
//    dS  = (S_max - S_min) / N
//
//  S_min defaults to 0 for standard options. For down-barrier options the
//  lower boundary is the barrier H > 0, so pass S_min = H explicitly:
//    Grid grid(S_max, T, N, M, H);   // domain [H, S_max]
//
//  Only one time slice is stored at a time (backward solve overwrites it).
//  After solve(), V_[i] holds V(t=0, S_i) for all i.
// =============================================================================

class Grid {
public:
    // S_min defaults to 0 — all existing Grid(S_max, T, N, M) calls still work.
    Grid(double S_max, double T, int N, int M, double S_min = 0.0);

    // -- Grid parameters -------------------------------------------------------
    int    N()     const { return N_; }       // number of spatial steps
    int    M()     const { return M_; }       // number of time steps
    double dS()    const { return dS_; }      // spatial step size
    double dt()    const { return dt_; }      // time step size
    double T()     const { return T_; }       // maturity
    double S_min() const { return S_min_; }   // lower spatial boundary
    double S_max() const { return S_max_; }   // upper spatial boundary

    double S(int i) const { return S_min_ + i * dS_; }   // spot price at node i
    double t(int m) const { return m * dt_; }             // time at step m

    // -- Solution access -------------------------------------------------------
    double& operator[](int i)       { return V_[i]; }
    double  operator[](int i) const { return V_[i]; }

    const std::vector<double>& values() const { return V_; }
          std::vector<double>& values()       { return V_; }

    // -- Initialisation --------------------------------------------------------

    // Set V_[i] = payoff(S_i) for all i  ->  terminal condition at t = T.
    void initialize(const Option& opt);

    // Enforce V_[0] and V_[N] from the option's boundary conditions at t_current.
    // Called by the solver at each time step before assembling the linear system.
    void apply_bc(const Option& opt, double t_current, double r);

private:
    double S_min_, S_max_, T_;
    int    N_, M_;
    double dS_, dt_;
    std::vector<double> V_;   // current time slice: V_[i] = V(t_current, S_i)
};
