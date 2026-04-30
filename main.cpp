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
//  Display helpers — European options (with BS analytical comparison)
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

// Run a solver on a European option, compute all Greeks, and print results.
static void run_option(Solver& solver, const Option& opt,
                       const BSParams& params,
                       const std::vector<double>& spots) {
    Grid grid(S_MAX, T, N, M);
    solver.solve(grid, opt, params);

    GreeksCalculator g_calc(solver);
    Greeks g = g_calc.compute(grid, opt, params);

    std::cout << "  [" << solver.name() << "]\n\n";
    for (double S : spots)
        print_spot(grid, g, params, S, opt.strike(), opt.is_call());
}

// =============================================================================
//  Display helpers — Exotic options (FD price only, no BS reference)
// =============================================================================

// Solve and print FD prices at the given spots.
// S_max_ allows passing a custom upper boundary (e.g., barrier = H instead of 3*K).
static void run_exotic(Solver& solver, const Option& opt,
                       const BSParams& params, double S_max_,
                       const std::vector<double>& spots) {
    Grid grid(S_max_, T, N, M);
    solver.solve(grid, opt, params);

    std::cout << "  [" << solver.name() << "]\n\n";
    std::cout << std::fixed << std::setprecision(5);
    for (double S0 : spots) {
        int i0 = static_cast<int>(std::round(S0 / grid.dS()));
        i0 = std::max(1, std::min(i0, grid.N() - 1));
        std::cout << "  S = " << std::setw(8) << S0
                  << "   V_FD = " << std::setw(9) << grid[i0] << "\n";
    }
    std::cout << "\n";
}

// Solve European and American puts and print a side-by-side comparison showing
// the early exercise premium = American - European and the intrinsic value.
//
// Expected observations:
//   - American >= European at every S (the extra right is never negative).
//   - American >= max(K - S, 0) at every S (the option is never below intrinsic).
//   - Deep ITM (S << K): large early exercise premium — the holder should exercise.
//   - OTM (S >> K): zero premium — early exercise has no value.
static void run_american_vs_european(Solver& solver, const BSParams& params,
                                     const std::vector<double>& spots) {
    EuropeanPut eu(K);
    AmericanPut am(K);

    Grid grid_eu(S_MAX, T, N, M);
    Grid grid_am(S_MAX, T, N, M);
    solver.solve(grid_eu, eu, params);
    solver.solve(grid_am, am, params);

    std::cout << "  [" << solver.name() << "] — European vs American\n\n";
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "  " << std::setw(8)  << "S"
              << std::setw(13) << "European"
              << std::setw(13) << "American"
              << std::setw(13) << "EE Premium"
              << std::setw(13) << "Intrinsic"
              << "\n";
    std::cout << "  " << std::string(58, '-') << "\n";
    for (double S0 : spots) {
        int i0 = static_cast<int>(std::round(S0 / grid_eu.dS()));
        i0 = std::max(1, std::min(i0, N - 1));
        const double S_node = i0 * grid_eu.dS();   // actual grid node (may differ from S0)
        const double vEu    = grid_eu[i0];
        const double vAm    = grid_am[i0];
        const double intr   = std::max(K - S_node, 0.0);   // intrinsic at the grid node
        std::cout << "  " << std::setw(8)  << S_node
                  << std::setw(13) << vEu
                  << std::setw(13) << vAm
                  << std::setw(13) << (vAm - vEu)
                  << std::setw(13) << intr
                  << "\n";
    }
    std::cout << "\n";
}

// Look up the grid value at the node closest to S0.
// Works for any grid (uniform or shifted) since it uses grid.dS() and grid.S_min().
static double lookup(const Grid& g, double S0) {
    int i = static_cast<int>(std::round((S0 - g.S_min()) / g.dS()));
    i = std::max(1, std::min(i, g.N() - 1));
    return g[i];
}

// Display KO price, KI price (= Vanilla - KO), and parity check (should = Vanilla).
//
// Expected observations:
//   KO <= Vanilla  : the knock-out feature can only hurt.
//   KI >= 0        : the knock-in feature can only help (relative to nothing).
//   KO + KI        : must equal Vanilla exactly (parity, no model assumptions).
//   Discount near H: largest when S is close to H (high probability of touching).
static void run_barrier(Solver& solver,
                        const Option& vanilla_opt, Grid& grid_van,
                        const Option& ko_opt,      Grid& grid_ko,
                        const std::vector<double>& spots,
                        const std::string& label) {
    solver.solve(grid_van, vanilla_opt, {r, sigma, T});
    solver.solve(grid_ko,  ko_opt,      {r, sigma, T});

    std::cout << "  [" << solver.name() << "]  " << label << "\n\n";
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "  " << std::setw(9)  << "S"
              << std::setw(12) << "Vanilla"
              << std::setw(12) << "KO"
              << std::setw(12) << "KI=Van-KO"
              << std::setw(12) << "KO+KI"
              << "\n";
    std::cout << "  " << std::string(57, '-') << "\n";
    for (double S0 : spots) {
        const double vVan = lookup(grid_van, S0);
        const double vKO  = lookup(grid_ko,  S0);
        const double vKI  = vVan - vKO;
        std::cout << "  " << std::setw(9)  << S0
                  << std::setw(12) << vVan
                  << std::setw(12) << vKO
                  << std::setw(12) << vKI
                  << std::setw(12) << (vKO + vKI)
                  << "\n";
    }
    std::cout << "\n";
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

    // =========================================================================
    //  European options — validated against closed-form Black-Scholes
    // =========================================================================

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

    // =========================================================================
    //  Exotic options — FD prices only (no closed-form reference)
    // =========================================================================

    // -- Digital Call ----------------------------------------------------------
    //
    //  dS = S_MAX / N is passed to the constructor so the sigmoid width adapts
    //  automatically to the grid resolution.
    //  Parity check: DigitalCall(S) + DigitalPut(S) should equal e^{-rT} ~ 0.951.
    {
        const double dS = S_MAX / N;
        DigitalCall dig(K, dS);
        const std::vector<double> dig_spots = {70.0, 95.0, 100.0, 105.0, 130.0};
        std::cout << "== Digital Call  (K=" << K << ", T=" << T
                  << ", sigma=" << sigma << ", r=" << r << ") ==\n\n";
        run_exotic(cn,  dig, params, S_MAX, dig_spots);
        run_exotic(rcn, dig, params, S_MAX, dig_spots);
    }

    // -- Digital Put -----------------------------------------------------------
    {
        const double dS = S_MAX / N;
        DigitalPut dig(K, dS);
        const std::vector<double> dig_spots = {130.0, 105.0, 100.0, 95.0, 70.0};
        std::cout << "== Digital Put   (K=" << K << ", T=" << T
                  << ", sigma=" << sigma << ", r=" << r << ") ==\n\n";
        run_exotic(cn,  dig, params, S_MAX, dig_spots);
        run_exotic(rcn, dig, params, S_MAX, dig_spots);
    }

    // =========================================================================
    //  Barrier options — all 8 types (4 KO solved by PDE, 4 KI via parity)
    //  KO + KI = Vanilla  →  last column "KO+KI" must equal "Vanilla"
    // =========================================================================

    // -- Up barriers  (H = 150 > K = 100,  S_max = H) -------------------------
    {
        constexpr double H_UP = 150.0;
        const std::vector<double> spots = {70.0, 100.0, 130.0};

        std::cout << "== Up-and-Out Call  (K=" << K << ", H=" << H_UP << ") ==\n\n";
        for (Solver* s : {(Solver*)&cn, (Solver*)&rcn}) {
            EuropeanCall  van(K);
            UpAndOutCall  ko(K, H_UP);
            Grid gv(S_MAX, T, N, M);
            Grid gk(H_UP,  T, N, M);   // S_max = H
            run_barrier(*s, van, gv, ko, gk, spots, "UOC");
        }

        std::cout << "== Up-and-Out Put   (K=" << K << ", H=" << H_UP << ") ==\n\n";
        for (Solver* s : {(Solver*)&cn, (Solver*)&rcn}) {
            EuropeanPut  van(K);
            UpAndOutPut  ko(K, H_UP);
            Grid gv(S_MAX, T, N, M);
            Grid gk(H_UP,  T, N, M);   // S_max = H
            run_barrier(*s, van, gv, ko, gk, spots, "UOP");
        }
    }

    // -- Down barriers  (H = 80 < K = 100,  S_min = H) ------------------------
    {
        constexpr double H_DN = 80.0;
        const std::vector<double> spots = {85.0, 100.0, 120.0};

        std::cout << "== Down-and-Out Call  (K=" << K << ", H=" << H_DN << ") ==\n\n";
        for (Solver* s : {(Solver*)&cn, (Solver*)&rcn}) {
            EuropeanCall   van(K);
            DownAndOutCall ko(K, H_DN);
            Grid gv(S_MAX, T, N, M);
            Grid gk(S_MAX, T, N, M, H_DN);   // S_min = H
            run_barrier(*s, van, gv, ko, gk, spots, "DOC");
        }

        std::cout << "== Down-and-Out Put   (K=" << K << ", H=" << H_DN << ") ==\n\n";
        for (Solver* s : {(Solver*)&cn, (Solver*)&rcn}) {
            EuropeanPut   van(K);
            DownAndOutPut ko(K, H_DN);
            Grid gv(S_MAX, T, N, M);
            Grid gk(S_MAX, T, N, M, H_DN);   // S_min = H
            run_barrier(*s, van, gv, ko, gk, spots, "DOP");
        }
    }

    // -- American Put ----------------------------------------------------------
    //
    //  Expected: American >= European; early exercise premium is largest for
    //  deep ITM puts (S = 70) where exercising now and earning r*K dominates.
    //  For OTM puts (S = 130) the premium should be near zero.
    {
        const std::vector<double> am_spots = {50.0, 70.0, 85.0, 100.0, 115.0, 130.0};
        std::cout << "== American Put  (K=" << K << ", T=" << T
                  << ", sigma=" << sigma << ", r=" << r << ") ==\n\n";
        run_american_vs_european(cn,  params, am_spots);
        run_american_vs_european(rcn, params, am_spots);
    }

    return 0;
}
