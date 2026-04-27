#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

#include "Option.hpp"
#include "Grid.hpp"
#include "Solver.hpp"
#include "CrankNicolson.hpp"
#include "ReducedCN.hpp"
#include "Greeks.hpp"
#include "BSFormula.hpp"

// =============================================================================
//  Problem parameters
// =============================================================================
static constexpr double S_MAX = 300.0;   // upper spatial boundary (= 3 * K)
static constexpr double T     = 1.0;     // maturity in years
static constexpr double r     = 0.05;    // risk-free rate (annualised)
static constexpr double sigma = 0.20;    // volatility (annualised)
static constexpr double K     = 100.0;   // strike price
static constexpr int    N     = 1000;    // number of spatial steps
static constexpr int    M     = 1000;    // number of time steps

// =============================================================================
//  Display helpers
// =============================================================================

// Classify a spot as ITM / ATM / OTM.
// A 2% band around K is used as the ATM region to avoid false classification
// due to the discrete grid spacing.
static std::string moneyness(double S, double strike, bool is_call) {
    const bool itm = is_call ? (S > strike * 1.02) : (S < strike * 0.98);
    const bool otm = is_call ? (S < strike * 0.98) : (S > strike * 1.02);
    if (itm) return "[ITM]";
    if (otm) return "[OTM]";
    return "[ATM]";
}

// Print FD price and all Greeks for a given spot, side-by-side with the
// closed-form BS values so the approximation error is immediately visible.
static void print_spot(const Grid& grid, const Greeks& g,
                       const BSParams& p, double S0, double strike, bool is_call) {
    // Find the closest interior grid node to S0.
    int i0 = static_cast<int>(std::round(S0 / grid.dS()));
    i0 = std::max(1, std::min(i0, grid.N() - 1));

    const double tau  = p.T;   // evaluated at t=0, so tau = T
    const double V_fd = grid[i0];
    const double V_bs = is_call ? BSAnalytical::call_price(S0, strike, p.r, p.sigma, tau)
                                : BSAnalytical::put_price (S0, strike, p.r, p.sigma, tau);

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "  S = " << std::setw(6) << S0 << "  " << moneyness(S0, strike, is_call) << "\n";
    std::cout << "    Price  FD = " << std::setw(9) << V_fd
              << "   BS = " << std::setw(9) << V_bs
              << "   err = " << std::setw(9) << std::abs(V_fd - V_bs) << "\n";
    std::cout << "    Delta  FD = " << std::setw(9) << g.delta[i0]
              << "   BS = " << std::setw(9)
              << (is_call ? BSAnalytical::call_delta(S0, strike, p.r, p.sigma, tau)
                          : BSAnalytical::put_delta (S0, strike, p.r, p.sigma, tau)) << "\n";
    std::cout << "    Gamma  FD = " << std::setw(9) << g.gamma[i0]
              << "   BS = " << std::setw(9)
              << BSAnalytical::gamma(S0, strike, p.r, p.sigma, tau) << "\n";
    std::cout << "    Theta  FD = " << std::setw(9) << g.theta[i0] / 365.0
              << "   BS = " << std::setw(9)
              << (is_call ? BSAnalytical::call_theta(S0, strike, p.r, p.sigma, tau) / 365.0
                          : BSAnalytical::put_theta (S0, strike, p.r, p.sigma, tau) / 365.0)
              << "  (daily)\n";
    std::cout << "    Vega   FD = " << std::setw(9) << g.vega[i0]
              << "   BS = " << std::setw(9)
              << BSAnalytical::vega(S0, strike, p.r, p.sigma, tau) << "\n";
    std::cout << "    Rho    FD = " << std::setw(9) << g.rho[i0]
              << "   BS = " << std::setw(9)
              << (is_call ? BSAnalytical::call_rho(S0, strike, p.r, p.sigma, tau)
                          : BSAnalytical::put_rho (S0, strike, p.r, p.sigma, tau)) << "\n";
    std::cout << "\n";
}

// Run a solver on an option, compute all Greeks, and print results for each spot.
static void run_option(Solver& solver, const Option& opt,
                       const BSParams& params,
                       const std::vector<double>& spots) {
    Grid grid(S_MAX, T, N, M);
    solver.solve(grid, opt, params);

    // GreeksCalculator holds a reference to solver for bump-and-reprice Greeks.
    GreeksCalculator g_calc(solver);
    Greeks g = g_calc.compute(grid, opt, params);

    std::cout << "  [" << solver.name() << "]\n\n";
    for (double S : spots)
        print_spot(grid, g, params, S, opt.strike(), opt.is_call());
}

// =============================================================================
//  Convergence analysis (uncomment to run)
// =============================================================================
//
//  Sweep N in {50, 100, 200, 500, 1000} with M = N and compute the L-inf
//  error against BSAnalytical::call_price over all interior nodes.
//  Plot log2(error) vs log2(N): expected slope ≈ -2 (CN is 2nd-order).
//  Compare CrankNicolson (uniform S-grid) vs ReducedCN (log-uniform grid).
//
//  static double max_error(Solver& solver, int Ni, const BSParams& p) {
//      EuropeanCall call(K);
//      Grid grid(S_MAX, T, Ni, Ni);
//      solver.solve(grid, call, p);
//      double err = 0.0;
//      for (int i = 1; i < Ni; ++i) {
//          double ref = BSAnalytical::call_price(grid.S(i), K, p.r, p.sigma, p.T);
//          err = std::max(err, std::abs(grid[i] - ref));
//      }
//      return err;
//  }
//
//  CrankNicolson cn;  ReducedCN rcn;
//  for (int Ni : {50, 100, 200, 500, 1000})
//      std::cout << "N=" << Ni
//                << "  CN="  << max_error(cn,  Ni, params)
//                << "  RCN=" << max_error(rcn, Ni, params) << "\n";

// =============================================================================
//  main
// =============================================================================
int main() {
    const BSParams params{r, sigma, T};

    CrankNicolson cn;    // standard CN on uniform S-grid
    ReducedCN     rcn;   // CN on log-uniform grid (constant coefficients)

    // Representative spots: OTM, ATM, ITM
    const std::vector<double> call_spots = { 70.0, 100.0, 130.0};
    const std::vector<double> put_spots  = {130.0, 100.0,  70.0};

    // -- European Call ---------------------------------------------------------
    {
        EuropeanCall call(K);
        std::cout << "== European Call  (K=" << K << ", T=" << T
                  << ", sigma=" << sigma << ", r=" << r << ") ==\n\n";
        run_option(cn,  call, params, call_spots);
        run_option(rcn, call, params, call_spots);
    }

    // -- European Put ----------------------------------------------------------
    {
        EuropeanPut put(K);
        std::cout << "== European Put   (K=" << K << ", T=" << T
                  << ", sigma=" << sigma << ", r=" << r << ") ==\n\n";
        run_option(cn,  put, params, put_spots);
        run_option(rcn, put, params, put_spots);
    }

    return 0;
}
