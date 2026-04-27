# Black-Scholes PDE Pricer — C++

![Status](https://img.shields.io/badge/status-Phase%200%20complete-success)

A C++ implementation of European option pricing by solving the Black-Scholes PDE numerically,
with two independent finite-difference schemes, full Greeks computation, and validation against
closed-form analytical formulas.

Initial commit built in collaboration with [Vithuson Vaithilingam](https://github.com/VithuVa).

---

## The problem

The Black-Scholes PDE governs the fair value $V(t, S)$ of a European option:

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S} - rV = 0$$

with terminal condition $V(T, S) = \text{Payoff}(S)$ and boundary conditions depending on the
option type. The PDE is solved **backwards in time** on a discrete $(S, t)$ grid, from maturity
$t = T$ down to today $t = 0$.

---

## Numerical methods

### Crank-Nicolson on a uniform S-grid

Standard centred finite differences on $S_i = i \cdot \Delta S$ yield a spatial operator:

$$\mathcal{L}[V]_i = a_i V_{i-1} + b_i V_i + c_i V_{i+1}$$

where $a_i, b_i, c_i$ depend on $S_i$ (variable coefficients). The Crank-Nicolson scheme
averages the implicit and explicit levels:

$$\frac{V_i^{n-1} - V_i^n}{\Delta t} = \frac{1}{2}\left(\mathcal{L}[V^{n-1}]_i + \mathcal{L}[V^n]_i\right)$$

This produces a **tridiagonal linear system** per time step, solved in $O(N)$ via the Thomas
algorithm. Properties: second-order in both time $O(\Delta t^2)$ and space $O(\Delta S^2)$,
unconditionally stable (von Neumann analysis).

### ReducedCN — Crank-Nicolson in log-price space

The substitution $x = \ln S$ transforms the BS PDE into a **constant-coefficient** PDE:

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 \frac{\partial^2 V}{\partial x^2} + \kappa \frac{\partial V}{\partial x} - rV = 0, \qquad \kappa = r - \tfrac{1}{2}\sigma^2$$

Key advantages over the S-space scheme:
- Coefficients $a, b, c$ are **constant** across all nodes — no growth for large $S$
- The log-uniform grid $\Delta x = \text{const}$ gives finer relative resolution near the strike
- $x = \ln S$ is the natural variable of the GBM dynamics (Itô's lemma)

After solving, $V(x_i)$ is linearly interpolated back onto the uniform $S$-grid for output.

---

## Greeks

| Greek | Method | Cost |
|-------|--------|------|
| $\Delta = \partial V / \partial S$ | Central FD on solved grid | Free |
| $\Gamma = \partial^2 V / \partial S^2$ | Central FD on solved grid | Free |
| $\Theta = \partial V / \partial t$ | PDE identity: $\Theta = -\mathcal{L}[V]$ | Free |
| $\nu = \partial V / \partial \sigma$ | Bump-and-reprice $\sigma \pm \varepsilon$ | 2 solver calls |
| $\rho = \partial V / \partial r$ | Bump-and-reprice $r \pm \varepsilon$ | 2 solver calls |

Theta is recovered for free from the PDE itself: since $\partial V / \partial t = -\mathcal{L}[V]$,
Theta is just the spatial operator evaluated at the solution — no time differencing needed.

All Greeks are validated against closed-form Black-Scholes formulas (see `BSFormula.hpp`).

---

## Project structure

```
BS-Pricer/
├── main.cpp                # Driver: runs both solvers, prints FD vs analytical comparison
│
├── Option.hpp              # Abstract Option interface + EuropeanCall, EuropeanPut
│
├── Grid.hpp / Grid.cpp     # Uniform (S, t) grid — storage, initialisation, boundary conditions
│
├── Solver.hpp              # Abstract Solver interface + BSParams struct
│
├── CrankNicolson.hpp/.cpp  # CN scheme on uniform S-grid (variable coefficients)
├── ReducedCN.hpp/.cpp      # CN scheme on log-uniform x-grid (constant coefficients)
│
├── Greeks.hpp / Greeks.cpp # Delta, Gamma, Theta (FD spatial) + Vega, Rho (bump-and-reprice)
├── BSFormula.hpp/.cpp      # Closed-form BS pricing and Greeks (analytical reference)
│
├── Thomas.hpp / Thomas.cpp # O(N) tridiagonal solver (Thomas algorithm)
│
└── Makefile
```

---

## Build

Requires **GCC** with C++17:

```bash
make
./pricer
```

Or directly:

```bash
g++ -std=c++17 -O2 Grid.cpp CrankNicolson.cpp ReducedCN.cpp Greeks.cpp BSFormula.cpp Thomas.cpp main.cpp -o pricer
./pricer
```

---

## Sample output

```
== European Call  (K=100, T=1, sigma=0.2, r=0.05) ==

  [Crank-Nicolson]

  S = 70.00000  [OTM]
    Price  FD =   0.43398   BS =   0.44145   err =   0.00747
    Delta  FD =   0.07487   BS =   0.07588
    Gamma  FD =   0.01011   BS =   0.01020
    Theta  FD =  -0.00336   BS =  -0.00341  (daily)
    Vega   FD =   9.87998   BS =   9.99690
    Rho    FD =   4.79958   BS =   4.86983

  S = 100.00000  [ATM]
    Price  FD =  10.38696   BS =  10.45058   err =   0.06362
    ...
```

Errors are of order $10^{-2}$ for $N = M = 1000$, consistent with the $O(\Delta S^2, \Delta t^2)$
convergence of the Crank-Nicolson scheme.

---

## Key concepts

**Thomas algorithm** — solves a tridiagonal system $Ax = d$ in $O(N)$ by forward elimination
and back substitution, replacing the $O(N^3)$ cost of general Gaussian elimination.

**Crank-Nicolson** — implicit-explicit averaging achieving $O(\Delta t^2, \Delta S^2)$ accuracy
while remaining unconditionally stable, unlike the explicit scheme which requires
$\Delta t = O(\Delta S^2)$.

**Log-price change of variable** — $x = \ln S$ maps variable-coefficient BS into a
constant-coefficient PDE, eliminating the $S^2$ growth of the diffusion term and producing
a geometrically-spaced grid with better resolution near the strike.

**PDE identity for Theta** — rather than finite-differencing in time, Theta is read directly
from the spatial operator $\mathcal{L}[V]$ via the BS PDE itself, which is both cheaper
(no extra time slice) and more accurate.

**Polymorphic design** — `Option`, `Solver`, and `GreeksCalculator` are decoupled via abstract
interfaces. Swapping `CrankNicolson` for `ReducedCN` (or any future scheme) requires no
changes to the Greeks or output code.

---

## Roadmap

This is **phase 0** — European options, two FD schemes, full Greeks, analytical validation with a total new architecture.

Next steps:
- **Exotic options** — digital (binary) payoffs, knock-out barrier options, American puts (PSOR solver)
- **SDL2 visualisation** — interactive window displaying $V(S)$ and Greek curves in real time, with keyboard controls to adjust $\sigma$, $r$, $T$ on the fly

