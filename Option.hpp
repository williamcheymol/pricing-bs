#pragma once
#include <algorithm>
#include <cmath>

// =============================================================================
//  Option.hpp  —  Abstract option interface + European Call / Put
// =============================================================================
//
//  Time convention : t in [0, T]
//    t = 0  ->  today          (value to compute)
//    t = T  ->  maturity       (terminal condition = payoff)
//
//  The time-to-maturity is  tau = T - t  (decreases as t increases).
//
//  Every concrete option must implement three boundary conditions that
//  constrain the PDE solver on the spatial domain S in [0, S_max] :
//    payoff  : terminal condition  V(T, S)
//    bc_lower: left  boundary      V(t, 0)
//    bc_upper: right boundary      V(t, S_max)
// =============================================================================

class Option {
public:
    virtual ~Option() = default;

    // Terminal condition: V(T, S) = payoff at expiry
    virtual double payoff(double S) const = 0;

    // Left boundary: V(t, S=0)
    // Called at every time step during the backward solve.
    virtual double bc_lower(double t, double T, double r) const = 0;

    // Right boundary: V(t, S=S_max)
    // S_max is passed explicitly — required for the call (see bug note below).
    virtual double bc_upper(double S_max, double t, double T, double r) const = 0;

    virtual const char* name()    const = 0;
    virtual double      strike()  const = 0;
    virtual bool        is_call() const = 0;
};

// =============================================================================
//  European Call  :  V(T, S) = max(S - K, 0)
// =============================================================================
class EuropeanCall : public Option {
public:
    explicit EuropeanCall(double K) : K_(K) {}

    // Payoff: the holder profits from every dollar above K.
    double payoff(double S) const override {
        return std::max(S - K_, 0.0);
    }

    // S -> 0 : the underlying can no longer rise above K, so the call is worthless.
    double bc_lower(double /*t*/, double /*T*/, double /*r*/) const override {
        return 0.0;
    }

    // S -> S_max : deeply ITM call behaves like a forward (Delta ≈ 1).
    //   C(S_max, t) ≈ S_max - K * e^{-r*tau}   where tau = T - t
    //
    // Common bug: returning K*exp(r*(T-t)) (wrong sign + missing S_max).
    double bc_upper(double S_max, double t, double T, double r) const override {
        return S_max - K_ * std::exp(-r * (T - t));
    }

    const char* name()    const override { return "European Call"; }
    double      strike()  const override { return K_; }
    bool        is_call() const override { return true; }

private:
    double K_;
};

// =============================================================================
//  European Put  :  V(T, S) = max(K - S, 0)
// =============================================================================
class EuropeanPut : public Option {
public:
    explicit EuropeanPut(double K) : K_(K) {}

    // Payoff: the holder profits from every dollar below K.
    double payoff(double S) const override {
        return std::max(K_ - S, 0.0);
    }

    // S -> 0 : exercise is certain, so the put equals the discounted strike.
    //   P(t, 0) = K * e^{-r*tau}
    double bc_lower(double t, double T, double r) const override {
        return K_ * std::exp(-r * (T - t));
    }

    // S -> S_max : deeply OTM put, essentially worthless.
    double bc_upper(double /*S_max*/, double /*t*/, double /*T*/, double /*r*/) const override {
        return 0.0;
    }

    const char* name()    const override { return "European Put"; }
    double      strike()  const override { return K_; }
    bool        is_call() const override { return false; }

private:
    double K_;
};
