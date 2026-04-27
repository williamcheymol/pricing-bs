#pragma once
#include "Solver.hpp"
#include "Thomas.hpp"
#include <vector>

// =============================================================================
//  ReducedCN.hpp  —  Crank-Nicolson on the log-price reduced PDE
// =============================================================================
//
//  Change of variable:  x = ln(S),  S = e^x
//
//  Starting from the standard BS PDE in S-space:
//    dV/dt  +  1/2*sigma^2*S^2 * V_SS  +  r*S * V_S  -  r*V  =  0
//
//  With  S*V_S = V_x  and  S^2*V_SS = V_xx - V_x  (chain rule + Ito), this
//  becomes the *reduced* (constant-coefficient) PDE in x-space:
//
//    dV/dt  +  1/2*sigma^2 * V_xx  +  kappa * V_x  -  r*V  =  0
//
//  where:
//    kappa = r - 1/2*sigma^2   (risk-neutral log-drift, from Ito's lemma)
//
//  -- Advantages over CrankNicolson (S-space) --------------------------------
//
//    1. CONSTANT coefficients: a, b, c do not depend on S_i.
//       -> uniform tridiagonal matrix, no large coefficients for large S.
//
//    2. Log-uniform grid: dx = const  <=>  dS/S = const (relative spacing).
//       -> finer resolution near the strike (ATM region), where V changes most.
//
//    3. Natural variable: x = ln(S) is exactly what drives a GBM
//       (dS/S = mu*dt + sigma*dW  <=>  dx = (mu - sigma^2/2)*dt + sigma*dW).
//
//  -- Output interpolation ---------------------------------------------------
//
//    The solve runs on the log-grid {x_i = ln(dS) + i*dx}.
//    At the end, V(x_i) is linearly interpolated onto the uniform S-grid
//    of the Grid object so the result is compatible with GreeksCalculator.
//
//    Note: this interpolation slightly degrades Gamma and Theta accuracy
//    (second derivatives amplify the interpolation error), while Price and
//    Delta remain as accurate as the standard CN scheme.
//
// =============================================================================

class ReducedCN : public Solver {
public:
    void        solve(Grid& grid, const Option& opt, const BSParams& p) override;
    const char* name() const override { return "ReducedCN (log-space)"; }

private:
    // Advance one time step on the log-grid: V^n -> V^{n-1}.
    void step(std::vector<double>& V,
              const Option&        opt,
              const BSParams&      p,
              double               x_min,
              double               dx,
              int                  Nx,
              double               dt,
              double               S_max,
              double               t_next);
};
