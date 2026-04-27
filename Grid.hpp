#pragma once
#include <vector>
#include "Option.hpp"

// =============================================================================
//  Grid.hpp  —  Uniform spatial / temporal discretisation grid
// =============================================================================
//
//  Uniform grid over (S, t) :
//    S_i = i * dS,   i = 0, 1, ..., N     (N+1 spatial nodes)
//    t_m = m * dt,   m = 0, 1, ..., M     (t=0 present, t=T maturity)
//
//  Only one time slice is stored at a time (backward solve overwrites it).
//  After solve(), V_[i] holds V(t=0, S_i) for all i.
//
//  Grid parameters are chosen so that S_max = 3*K is typically sufficient:
//  beyond 3K a call's Delta is ≈ 1 and the upper BC is tight analytically.
// =============================================================================

class Grid {
public:
    Grid(double S_max, double T, int N, int M);

    // -- Grid parameters -------------------------------------------------------
    int    N()     const { return N_; }       // number of spatial steps
    int    M()     const { return M_; }       // number of time steps
    double dS()    const { return dS_; }      // spatial step size
    double dt()    const { return dt_; }      // time step size
    double T()     const { return T_; }       // maturity
    double S_max() const { return S_max_; }   // upper spatial boundary

    double S(int i) const { return i * dS_; }   // spot price at node i
    double t(int m) const { return m * dt_; }   // time at step m

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
    double S_max_, T_;
    int    N_, M_;
    double dS_, dt_;
    std::vector<double> V_;   // current time slice: V_[i] = V(t_current, S_i)
};
