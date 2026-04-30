#include "Grid.hpp"

Grid::Grid(double S_max, double T, int N, int M, double S_min)
    : S_min_(S_min), S_max_(S_max), T_(T), N_(N), M_(M),
      dS_((S_max - S_min) / N),
      dt_(T / M),
      V_(N + 1, 0.0)
{}

// Fill the grid with the option payoff at maturity t = T.
// This is the terminal condition that bootstraps the backward solve.
void Grid::initialize(const Option& opt) {
    for (int i = 0; i <= N_; ++i)
        V_[i] = opt.payoff(S(i));
}

// Impose boundary conditions at the current time level t_current.
// The solver calls this before building each tridiagonal system so that
// V_[0] and V_[N_] hold the correct analytical BCs at the new time.
void Grid::apply_bc(const Option& opt, double t_current, double r) {
    V_[0]  = opt.bc_lower(t_current, T_, r);
    V_[N_] = opt.bc_upper(S_max_, t_current, T_, r);
}
