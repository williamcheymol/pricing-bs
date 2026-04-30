#pragma once
#include <algorithm>
#include <cmath>

// =============================================================================
//  Option.hpp  —  Abstract option interface + concrete payoff classes
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
//
//  Options implemented here:
//    European  : Call, Put
//    Digital   : Call, Put  (cash-or-nothing, discontinuous payoff)
//    Barrier   : Up-and-Out Call  (S_max must be set to H in the Grid)
//    American  : Put  (early exercise via the intrinsic() projection)
// =============================================================================

class Option {
public:
    virtual ~Option() = default;

    // Terminal condition: V(T, S) = payoff at expiry.
    virtual double payoff(double S) const = 0;

    // Left boundary: V(t, S=0).
    // Called at every time step during the backward solve.
    virtual double bc_lower(double t, double T, double r) const = 0;

    // Right boundary: V(t, S=S_max).
    // S_max is passed explicitly — required for the call (see EuropeanCall comment).
    virtual double bc_upper(double S_max, double t, double T, double r) const = 0;

    virtual const char* name()    const = 0;
    virtual double      strike()  const = 0;
    virtual bool        is_call() const = 0;

    // Intrinsic (immediate exercise) value at spot S.
    //
    // Used by both CN solvers to enforce the early-exercise constraint:
    //   V_i = max(V_i, intrinsic(S_i))   after each backward step.
    //
    // For European and path-dependent options: default returns 0.0, meaning
    // the constraint is never active and the projection costs nothing.
    // Only American options need to override this.
    virtual double intrinsic(double /*S*/) const { return 0.0; }
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

// =============================================================================
//  Digital Call (Cash-or-Nothing)  :  V(T, S) = 1_{S > K}
// =============================================================================
//
//  Economic intuition
//  ------------------
//  The digital call pays exactly 1 unit of cash if the underlying closes above
//  K — regardless of how far above. Unlike a vanilla call, the payout is flat:
//  once you are in the money, you have collected the full prize. The digital is
//  therefore a pure probability bet: under risk-neutral pricing,
//
//      Digital Call = e^{-rT} * P*(S_T > K) = e^{-rT} * N(d2).
//
//  Its Delta is the risk-neutral PDF of S_T, which peaks sharply near K — far
//  larger than a vanilla call's Delta near the strike.
//
//  Numerical challenge: Gibbs phenomenon
//  --------------------------------------
//  The payoff has a unit jump at S = K (step function). Centred finite
//  differences resolve smooth functions to O(dS^2) but produce ringing
//  oscillations near discontinuities — the Gibbs phenomenon. These ripples
//  contaminate Delta and Gamma near K and can cause negative option prices.
//  Remedies:
//    1. Payoff smoothing [implemented]: replace the step with a sigmoid of
//       width ~dS centred at K — eliminates the jump while keeping O(dS) error.
//    2. Rannacher time-stepping: use 2 fully implicit steps (theta=1) at the
//       start to damp high-frequency modes before switching to CN.
//    3. Non-uniform grid: cluster nodes near K to resolve the jump.
//
// =============================================================================
class DigitalCall : public Option {
public:
    // dS = spatial step of the grid — sets the sigmoid width for payoff smoothing.
    // Pass S_max / N from main.cpp so the smoothing adapts to the grid resolution.
    DigitalCall(double K, double dS) : K_(K), dS_(dS) {}

    // Smoothed payoff: sigmoid of width dS centred at K.
    //
    // Why not the raw step 1_{S > K}?
    // The step function has an infinite derivative at K. Centred FD approximates
    // smooth functions well but produces Gibbs oscillations near jumps that never
    // converge away. Replacing the step with a sigmoid of width ~dS eliminates
    // the discontinuity: the function is now C-infinity, the scheme converges
    // cleanly, and the pricing error introduced is O(dS) — same order as the
    // discretisation error already present.
    double payoff(double S) const override {
        return 1.0 / (1.0 + std::exp(-2.0 * (S - K_) / dS_));
    }

    // S -> 0 : certain loser, worthless.
    double bc_lower(double /*t*/, double /*T*/, double /*r*/) const override {
        return 0.0;
    }

    // S -> S_max : certain winner — discounted value of receiving 1 at maturity.
    //   V = e^{-r*tau}   (tau = T - t)
    // Note the minus sign: we DISCOUNT (bring future cash to present value).
    double bc_upper(double /*S_max*/, double t, double T, double r) const override {
        return std::exp(-r * (T - t));
    }

    const char* name()    const override { return "Digital Call"; }
    double      strike()  const override { return K_; }
    bool        is_call() const override { return true; }

private:
    double K_, dS_;
};

// =============================================================================
//  Digital Put (Cash-or-Nothing)  :  V(T, S) = 1_{S < K}
// =============================================================================
//
//  Economic intuition
//  ------------------
//  The digital put is the mirror of the digital call. It pays 1 if S_T < K.
//  Under put-call parity for digitals:
//
//      Digital Call + Digital Put = e^{-rT}    (one side always pays)
//
//  so Digital Put = e^{-rT} * N(-d2) = e^{-rT} * (1 - N(d2)).
//
//  Use this parity to check your implementation:
//  for any S, DigitalCall(S) + DigitalPut(S) should equal e^{-rT}.
//
// =============================================================================
class DigitalPut : public Option {
public:
    // dS = spatial step — same smoothing logic as DigitalCall.
    DigitalPut(double K, double dS) : K_(K), dS_(dS) {}

    // Smoothed payoff: 1 - sigmoid (mirror of DigitalCall).
    //
    // Parity check: DigitalCall::payoff(S) + DigitalPut::payoff(S) = 1 exactly,
    // which matches the discounted parity DC + DP = e^{-rT} at t = T (tau = 0).
    double payoff(double S) const override {
        return 1.0 - 1.0 / (1.0 + std::exp(-2.0 * (S - K_) / dS_));
    }

    // S -> 0 : certain winner — same discounted value as DigitalCall's bc_upper.
    double bc_lower(double t, double T, double r) const override {
        return std::exp(-r * (T - t));
    }

    // S -> S_max : certain loser, worthless.
    double bc_upper(double /*S_max*/, double /*t*/, double /*T*/, double /*r*/) const override {
        return 0.0;
    }

    const char* name()    const override { return "Digital Put"; }
    double      strike()  const override { return K_; }
    bool        is_call() const override { return false; }

private:
    double K_, dS_;
};

// =============================================================================
//  Barrier options — general principles
// =============================================================================
//
//  A barrier option is a vanilla option with an additional knock-out (KO) or
//  knock-in (KI) condition tied to a barrier level H:
//
//    Knock-Out : the option is CANCELLED if S ever touches H (pays 0 from then)
//    Knock-In  : the option only ACTIVATES if S ever touches H (else pays 0)
//
//  Fundamental parity (always true, no model assumptions):
//
//      Knock-Out  +  Knock-In  =  Vanilla
//
//  This means KI options are FREE once you have the KO prices:
//      KI(S) = Vanilla(S) - KO(S)
//  No extra PDE solve needed — computed in main.cpp by subtraction.
//
//  Four KO types, depending on barrier direction and option type:
//
//    Up-and-Out Call  (UOC) : H > K,  S_max = H,  bc_upper = 0
//    Up-and-Out Put   (UOP) : H > K,  S_max = H,  bc_upper = 0
//    Down-and-Out Call(DOC) : H < K,  S_min = H,  bc_lower = 0
//    Down-and-Out Put (DOP) : H < K,  S_min = H,  bc_lower = 0
//
//  Grid construction rules:
//    Up barriers   → Grid(H,     T, N, M)        S_max = H
//    Down barriers → Grid(S_max, T, N, M, H)     S_min = H
//
// =============================================================================

// =============================================================================
//  Up-and-Out Call (UOC)  :  max(S_T - K, 0)  if  max S_t < H,  else 0
// =============================================================================
//
//  Economic intuition: bullish view with a cap — you profit if S rises above K
//  but NOT all the way to H. Cheaper than a vanilla call because you give up
//  the payoff in strongly bullish scenarios. Discount = value of paths hitting H.
//
//  !! Grid must be constructed with S_max = H !!
// =============================================================================
class UpAndOutCall : public Option {
public:
    UpAndOutCall(double K, double H) : K_(K), H_(H) {}

    // Knocked out if S >= H, otherwise vanilla call payoff.
    double payoff(double S) const override {
        return (S >= H_) ? 0.0 : std::max(S - K_, 0.0);
    }

    // S -> 0 : call worthless regardless of barrier (same as vanilla call).
    double bc_lower(double /*t*/, double /*T*/, double /*r*/) const override {
        return 0.0;
    }

    // S = H : barrier touched — knocked out.
    double bc_upper(double /*S_max*/, double /*t*/, double /*T*/, double /*r*/) const override {
        return 0.0;
    }

    const char* name()    const override { return "Up-and-Out Call"; }
    double      strike()  const override { return K_; }
    bool        is_call() const override { return true; }
    double      barrier() const          { return H_; }

private:
    double K_, H_;
};

// =============================================================================
//  Up-and-Out Put (UOP)  :  max(K - S_T, 0)  if  max S_t < H,  else 0
// =============================================================================
//
//  Economic intuition: bearish view where you also bet the stock won't spike.
//  Cheaper than a vanilla put — you lose the payoff if S ever rallies to H.
//  Less common than UOC (a put buyer is usually bearish, so H >> S is less
//  of a constraint), but useful to hedge positions with an upside cap.
//
//  bc_lower subtlety: at S = 0, the underlying is stuck (GBM absorbs at 0)
//  so the up barrier H can NEVER be reached from there. The option behaves
//  exactly like a European put: bc_lower = K * e^{-r*tau}.
//
//  !! Grid must be constructed with S_max = H !!
// =============================================================================
class UpAndOutPut : public Option {
public:
    UpAndOutPut(double K, double H) : K_(K), H_(H) {}

    // Knocked out if S >= H, otherwise vanilla put payoff.
    double payoff(double S) const override {
        return (S >= H_) ? 0.0 : std::max(K_ - S, 0.0);
    }

    //       S -> 0, the up barrier can never be reached (GBM absorbs at 0).
    //       The option behaves like a European put — same bc_lower as EuropeanPut.
    double bc_lower(double t, double T, double r) const override {
        return K_ * std::exp(-r * (T - t));
    }

    // Barrier touched — Knock Out
    double bc_upper(double /*S_max*/, double /*t*/, double /*T*/, double /*r*/) const override {
        return 0.0;
    }

    const char* name()    const override { return "Up-and-Out Put"; }
    double      strike()  const override { return K_; }
    bool        is_call() const override { return false; }
    double      barrier() const          { return H_; }

private:
    double K_, H_;
};

// =============================================================================
//  Down-and-Out Call (DOC)  :  max(S_T - K, 0)  if  min S_t > H,  else 0
// =============================================================================
//
//  Economic intuition: you buy a call but lose everything if the stock falls
//  too far. Useful for hedging: a long stock + DOC protects the upside while
//  limiting the premium paid. If the stock crashes through H, both positions
//  suffer — but the call was cheap precisely because of that risk.
//
//  Typical setup: H < S_0 < K  (stock is between barrier and strike at inception).
//  The barrier is below the current spot — it only gets hit in a crash.
//
//  !! Grid must be constructed with S_min = H: Grid(S_max, T, N, M, H) !!
//  The lower boundary of the domain is H, not 0.
// =============================================================================
class DownAndOutCall : public Option {
public:
    DownAndOutCall(double K, double H) : K_(K), H_(H) {}

    // Knocked out if S <= H, otherwise vanilla call payoff.
    double payoff(double S) const override {
        return (S <= H_) ? 0.0 : std::max(S - K_, 0.0);
    }

    double bc_lower(double /*t*/, double /*T*/, double /*r*/) const override {
        return 0.0;
    }

    //       S -> S_max, deep ITM call far from barrier H << S_max.
    //       It behaves like a vanilla call — same bc_upper as EuropeanCall.
    double bc_upper(double S_max, double t, double T, double r) const override {
        return S_max - K_ * std::exp(-r * (T - t));
    }

    const char* name()    const override { return "Down-and-Out Call"; }
    double      strike()  const override { return K_; }
    bool        is_call() const override { return true; }
    double      barrier() const          { return H_; }

private:
    double K_, H_;
};

// =============================================================================
//  Down-and-Out Put (DOP)  :  max(K - S_T, 0)  if  min S_t > H,  else 0
// =============================================================================
//
//  Economic intuition: the most traded barrier option in practice. A bearish
//  investor buys a put but accepts that it's cancelled if the stock crashes all
//  the way to H — which seems contradictory, but H is typically set well below
//  the current spot. The logic: if the stock has already crashed to H, the put
//  is very deep ITM and the investor takes profits by closing the position early.
//  The DOP is significantly cheaper than a vanilla put.
//
//  Parity check: DOP + DIP = Vanilla Put
//  The DIP (Down-and-In Put) activates only on a crash — it's the "catastrophe
//  insurance" component. The DOP covers the normal bearish scenario.
//
//  !! Grid must be constructed with S_min = H: Grid(S_max, T, N, M, H) !!
// =============================================================================
class DownAndOutPut : public Option {
public:
    DownAndOutPut(double K, double H) : K_(K), H_(H) {}

    // Knocked out if S <= H, otherwise vanilla put payoff.
    double payoff(double S) const override {
        return (S <= H_) ? 0.0 : std::max(K_ - S, 0.0);
    }

    // bc_lower at S = H (the barrier).
    double bc_lower(double /*t*/, double /*T*/, double /*r*/) const override {
        return 0.0;
    }

    // S -> S_max, deep OTM put (same as European put upper BC).
    double bc_upper(double /*S_max*/, double /*t*/, double /*T*/, double /*r*/) const override {
        return 0.0;
    }

    const char* name()    const override { return "Down-and-Out Put"; }
    double      strike()  const override { return K_; }
    bool        is_call() const override { return false; }
    double      barrier() const          { return H_; }

private:
    double K_, H_;
};

// =============================================================================
//  American Put  :  V >= max(K - S, 0),  early exercise allowed at any t in [0, T]
// =============================================================================
//
//  Economic intuition
//  ------------------
//  A European put holder is forced to wait until T to collect max(K - S_T, 0).
//  An American put holder can exercise at any time: if the stock crashes deeply
//  below K, waiting is suboptimal because:
//    - The interest earned by receiving K now (and investing at rate r) may
//      exceed any further gain from the put staying alive.
//    - If S is already near 0, there is limited further downside to capture.
//
//  Early exercise is optimal in the region {S < S*(t)} where S*(t) is the
//  "optimal exercise boundary" — a curve that must be found simultaneously
//  with the option price (a free-boundary problem).
//
//  The American put is ALWAYS worth more than the European put:
//      American Put >= European Put    (extra right has non-negative value)
//      American Put >= max(K - S, 0)   (never worth less than immediate exercise)
//  The difference is the "early exercise premium."
//
//  For a call on a non-dividend-paying stock, early exercise is NEVER optimal:
//  holding the call always dominates (time value > 0). Hence,
//  American Call = European Call for non-dividend assets.
//
//  Projection method (linear complementarity)
//  -------------------------------------------
//  The CN scheme produces a candidate "continuation value" V*_i. The American
//  constraint V >= K - S is then enforced by a simple projection after each step:
//
//      V_i = max(V*_i,  K - S_i)
//
//  Nodes where V*_i < K - S_i lie in the exercise region (the PDE is inactive
//  there); nodes where V*_i >= K - S_i lie in the continuation region (the PDE
//  governs the price). This operator-splitting approach is exact for the CN
//  scheme and equivalent to PSOR with a single projection pass.
//
//  Possible extensions: extract the exercise boundary S*(t) by scanning for
//  the leftmost node where V_i > K - S_i; compare against the Barone-Adesi &
//  Whaley (1987) analytical approximation.
// =============================================================================
class AmericanPut : public Option {
public:
    explicit AmericanPut(double K) : K_(K) {}

    //       Same terminal condition as a European put — the payoff is collected
    //       at T only if the holder has not already exercised early.
    double payoff(double S) const override {
        return std::max(K_ - S, 0.0);
    }

    //       S -> 0, the put is deep ITM — early exercise is certain.
    //       The holder exercises immediately and collects K.
    //       Note: unlike the European put, the American put at S=0 is worth K,
    //       not K*e^{-r*tau} — there is no reason to wait for discounting.
    double bc_lower(double /*t*/, double /*T*/, double /*r*/) const override {
        return K_;
    }

    // S -> S_max >> K, deep OTM. Will the holder ever exercise?
    double bc_upper(double /*S_max*/, double /*t*/, double /*T*/, double /*r*/) const override {
        return 0.0;
    }

    //       intrinsic value — the payoff from exercising immediately at spot S.
    //       This is what the projection step compares against the continuation
    //       value to decide whether to exercise.
    double intrinsic(double S) const override {
        return K_ - S;
    }

    const char* name()    const override { return "American Put"; }
    double      strike()  const override { return K_; }
    bool        is_call() const override { return false; }

private:
    double K_;
};
