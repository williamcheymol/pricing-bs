#include "BSFormula.hpp"
#include <stdexcept>

// =============================================================================
//  BSFormula.cpp  —  Closed-form Black-Scholes pricing and Greeks
// =============================================================================

// -- Primitives ----------------------------------------------------------------

// Standard normal CDF via the complementary error function.
// Relation: N(x) = 0.5 * erfc(-x / sqrt(2))
// std::erfc is used instead of a custom approximation for maximum accuracy.
double BSAnalytical::N(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// d1 = [ ln(S/K) + (r + sigma^2/2) * tau ] / (sigma * sqrt(tau))
// Throws if tau <= 0: the formula is undefined at expiry (d1 -> +/-inf).
double BSAnalytical::d1(double S, double K, double r, double sigma, double tau) {
    if (tau <= 0.0)
        throw std::invalid_argument("tau must be > 0 (time to maturity)");
    return (std::log(S / K) + (r + 0.5 * sigma * sigma) * tau)
           / (sigma * std::sqrt(tau));
}

// d2 = d1 - sigma * sqrt(tau)
// The two quantities differ by exactly one standard deviation of the log-return.
double BSAnalytical::d2(double S, double K, double r, double sigma, double tau) {
    return d1(S, K, r, sigma, tau) - sigma * std::sqrt(tau);
}

// -- Pricing -------------------------------------------------------------------

// C = S * N(d1) - K * e^{-r*tau} * N(d2)
// Intuition: S*N(d1) is the risk-neutral expected payoff from receiving the
// stock; K*e^{-r*tau}*N(d2) is the discounted cost of paying the strike.
double BSAnalytical::call_price(double S, double K, double r, double sigma, double tau) {
    return S * N(d1(S, K, r, sigma, tau))
           - K * std::exp(-r * tau) * N(d2(S, K, r, sigma, tau));
}

// P = C - S + K * e^{-r*tau}   (put-call parity)
// Reuses call_price to avoid recomputing d1/d2.
double BSAnalytical::put_price(double S, double K, double r, double sigma, double tau) {
    return call_price(S, K, r, sigma, tau) - S + K * std::exp(-r * tau);
}

// -- Greeks --------------------------------------------------------------------

// Delta_call = N(d1)
// Probability (risk-neutral) that the call expires in the money.
double BSAnalytical::call_delta(double S, double K, double r, double sigma, double tau) {
    return N(d1(S, K, r, sigma, tau));
}

// Delta_put = N(d1) - 1
// Follows from put-call parity applied to deltas: Delta_put = Delta_call - 1.
double BSAnalytical::put_delta(double S, double K, double r, double sigma, double tau) {
    return N(d1(S, K, r, sigma, tau)) - 1.0;
}

// Gamma = phi(d1) / (S * sigma * sqrt(tau))
// Identical for calls and puts (put-call parity implies same convexity).
double BSAnalytical::gamma(double S, double K, double r, double sigma, double tau) {
    return phi(d1(S, K, r, sigma, tau)) / (S * sigma * std::sqrt(tau));
}

// Vega = S * phi(d1) * sqrt(tau)
// Identical for calls and puts.
// Measures the change in price for a 1-unit move in sigma (annualised).
double BSAnalytical::vega(double S, double K, double r, double sigma, double tau) {
    return S * phi(d1(S, K, r, sigma, tau)) * std::sqrt(tau);
}

// Theta_call = -S*phi(d1)*sigma/(2*sqrt(tau)) - r*K*e^{-r*tau}*N(d2)
// Both terms are negative for a long call: time decay always erodes option value.
double BSAnalytical::call_theta(double S, double K, double r, double sigma, double tau) {
    const double d1v = d1(S, K, r, sigma, tau);
    return - S * phi(d1v) * sigma / (2.0 * std::sqrt(tau))
           - r * K * std::exp(-r * tau) * N(d2(S, K, r, sigma, tau));
}

// Theta_put = -S*phi(d1)*sigma/(2*sqrt(tau)) + r*K*e^{-r*tau}*N(-d2)
// The second term is positive: deep ITM puts can have positive theta
// because the time value of money works in the holder's favour.
double BSAnalytical::put_theta(double S, double K, double r, double sigma, double tau) {
    const double d1v = d1(S, K, r, sigma, tau);
    return - S * phi(d1v) * sigma / (2.0 * std::sqrt(tau))
           + r * K * std::exp(-r * tau) * N(-d2(S, K, r, sigma, tau));
}

// Rho_call = +K * tau * e^{-r*tau} * N(d2)
// Always positive: a higher rate discounts the strike less, benefiting the call holder.
double BSAnalytical::call_rho(double S, double K, double r, double sigma, double tau) {
    return K * tau * std::exp(-r * tau) * N(d2(S, K, r, sigma, tau));
}

// Rho_put = -K * tau * e^{-r*tau} * N(-d2)
// Always negative: a higher rate discounts the strike received at exercise, hurting the put holder.
double BSAnalytical::put_rho(double S, double K, double r, double sigma, double tau) {
    return -K * tau * std::exp(-r * tau) * N(-d2(S, K, r, sigma, tau));
}
