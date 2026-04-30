#pragma once
#include <cmath>

// =============================================================================
//  BSFormula.hpp  —  Closed-form Black-Scholes pricing and Greeks
// =============================================================================
//
//  These analytical formulas serve as ground truth for validating the
//  finite-difference results.  They are exact for European options under
//  the constant-sigma BS model.
//
//  Standard notation:
//    S    : current spot price of the underlying
//    K    : strike price
//    tau  : time to maturity  (tau = T - t >= 0)
//    r    : risk-free rate (continuous, annualised)
//    sigma: implied volatility (annualised)
//
//  Intermediate quantities:
//    d1 = [ ln(S/K) + (r + sigma^2/2) * tau ] / (sigma * sqrt(tau))
//    d2 = d1 - sigma * sqrt(tau)
//
//  Pricing:
//    Call : C = S * N(d1) - K * e^{-r*tau} * N(d2)
//    Put  : P = C - S + K * e^{-r*tau}               (put-call parity)
//
//  Greeks (analytical derivatives):
//  +-------+-------------------------------------------------+
//  | Delta | call = N(d1)        put = N(d1) - 1             |
//  | Gamma | phi(d1) / (S*sigma*sqrt(tau))  (same call/put)  |
//  | Vega  | S * phi(d1) * sqrt(tau)        (same call/put)  |
//  | Theta | see formulas below (differ call vs put)         |
//  | Rho   | K*tau*e^{-r*tau} * N(d2) / N(-d2)              |
//  +-------+-------------------------------------------------+
// =============================================================================

struct BSAnalytical {

    // -- CDF and PDF of the standard normal ------------------------------------

    // N(x) = 0.5 * erfc(-x / sqrt(2))   (standard normal CDF)
    static double N(double x);

    // d1 and d2 as defined above.
    // Throws if tau <= 0 (undefined at expiry — handle the limit upstream).
    static double d1(double S, double K, double r, double sigma, double tau);
    static double d2(double S, double K, double r, double sigma, double tau);

    // -- Pricing ---------------------------------------------------------------

    // European call price: C = S*N(d1) - K*e^{-r*tau}*N(d2)
    static double call_price(double S, double K, double r, double sigma, double tau);

    // European put price via put-call parity: P = C - S + K*e^{-r*tau}
    // Reuses call_price — no need to recompute d1, d2.
    static double put_price(double S, double K, double r, double sigma, double tau);

    // -- Greeks ----------------------------------------------------------------

    // Delta: first derivative of V with respect to S.
    //   call: N(d1)         put: N(d1) - 1  (put-call parity on deltas)
    static double call_delta(double S, double K, double r, double sigma, double tau);
    static double put_delta (double S, double K, double r, double sigma, double tau);

    // Gamma: second derivative of V with respect to S. Identical for call and put.
    //   Gamma = phi(d1) / (S * sigma * sqrt(tau))
    static double gamma(double S, double K, double r, double sigma, double tau);

    // Vega: sensitivity to volatility. Identical for call and put.
    //   Vega = S * phi(d1) * sqrt(tau)
    //   Note: often quoted for a 1% move in vol (divide by 100).
    static double vega(double S, double K, double r, double sigma, double tau);

    // Theta: time decay (per year; divide by 365 for daily).
    //   call: -S*phi(d1)*sigma/(2*sqrt(tau)) - r*K*e^{-r*tau}*N(d2)
    //   put:  -S*phi(d1)*sigma/(2*sqrt(tau)) + r*K*e^{-r*tau}*N(-d2)
    static double call_theta(double S, double K, double r, double sigma, double tau);
    static double put_theta (double S, double K, double r, double sigma, double tau);

    // Rho: sensitivity to the risk-free rate.
    //   call: +K * tau * e^{-r*tau} * N(d2)
    //   put:  -K * tau * e^{-r*tau} * N(-d2)
    static double call_rho(double S, double K, double r, double sigma, double tau);
    static double put_rho (double S, double K, double r, double sigma, double tau);

    // -- Internal utility ------------------------------------------------------
    // Standard normal PDF: phi(x) = e^{-x^2/2} / sqrt(2*pi)
    static double phi(double x) {
        static constexpr double INV_SQRT2PI = 0.3989422804014327;
        return INV_SQRT2PI * std::exp(-0.5 * x * x);
    }
};
