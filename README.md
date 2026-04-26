# Black-Scholes PDE Pricing — C++

A C++ implementation of European option pricing by solving the Black-Scholes PDE numerically, with two finite-difference schemes and real-time SDL2 visualisation.

Initial commit built at **ENSIIE** (2024–2025), in collaboration with [Vithuson Vaithilingam](https://github.com/VithuVa).

---

## The problem

The Black-Scholes PDE governs the fair value $C(t, S)$ of a European option:

$$\frac{\partial C}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 C}{\partial S^2} + rS\frac{\partial C}{\partial S} - rC = 0$$

with terminal condition $C(T, S) = \text{Payoff}(S)$ and boundary conditions depending on the option type (call or put).

The PDE is solved backwards in time on a uniform $(S, t)$ grid using two independent finite-difference methods, then compared.

---

## Methods

### Full method — Crank-Nicolson

The Crank-Nicolson scheme averages the implicit and explicit operators at each time step:

$$\frac{C^{n+1}_i - C^n_i}{\Delta t} = \frac{1}{2}\left(\mathcal{L}C^{n+1}_i + \mathcal{L}C^n_i\right)$$

where $\mathcal{L}$ is the Black-Scholes spatial operator. This produces a **tridiagonal linear system** at each step, solved via the Thomas algorithm. Crank-Nicolson is second-order in both time and space and unconditionally stable.

### Reduced method — logarithmic change of variable

The substitution $x = \ln S$ transforms the Black-Scholes PDE into a **constant-coefficient heat equation**:

$$\frac{\partial u}{\partial \tau} = \frac{\partial^2 u}{\partial x^2}$$

on a uniform $x$-grid, solved with a fully implicit scheme. The solution is then mapped back to the original $(S, C)$ space by inverse transformation and interpolation. This formulation simplifies the spatial discretisation and eliminates the variable-coefficient stiffness of the direct approach.

---

## Project structure

```
bs-pricing/
├── main.cpp            # Entry point — runs both solvers, prints stats, launches SDL
│
├── Payoff.hpp          # Abstract Payoff base class + Call and Put
│
├── EDP.[cpp|hpp]       # PDE grid — discretisation, boundary conditions, solution storage
│
├── Res_full.[cpp|hpp]  # Crank-Nicolson solver (full PDE)
├── Res_red.[cpp|hpp]   # Implicit solver (reduced heat equation)
│
├── Thomas.[cpp|hpp]    # Thomas algorithm — O(n) tridiagonal solver
│
├── SDL.[cpp|hpp]       # SDL2 manager — initialisation, event loop, window ownership
├── window.[cpp|hpp]    # Window — renders curves pixel by pixel via SDL2
│
└── Doxyfile            # Doxygen configuration
```

---

## Build

Requires **GCC** and **SDL2**:

```bash
# Ubuntu/Debian
sudo apt install libsdl2-dev

# macOS
brew install sdl2
```

```bash
g++ -std=c++17 -O2 main.cpp EDP.cpp Res_full.cpp Res_red.cpp Thomas.cpp SDL.cpp window.cpp \
    $(sdl2-config --cflags --libs) -o prog
./prog
```

---

## Output

The program runs both solvers for **Call** and **Put** options and prints:

```
Prix au strike K (100) :
  C_full(0,K) = 5.174
  C_red(0,K)  = 5.171

--- Statistiques Erreur Absolue Call ---
  Min        : 0.000
  Max        : 0.032
  Moyenne    : 0.004
  Ecart-type : 0.006
```

Two SDL2 windows open per option type:
- **Window 1** — $C(0, S)$ curves: full method (blue) vs reduced (red)
- **Window 2** — Absolute error between the two methods (green)

---

## Key concepts

**Thomas algorithm** — solves a tridiagonal system $Ax = d$ in $O(n)$ by forward elimination and back substitution, avoiding the $O(n^3)$ cost of general Gaussian elimination.

**Crank-Nicolson** — implicit-explicit averaging that achieves $O(\Delta t^2, \Delta S^2)$ convergence while remaining unconditionally stable, unlike the explicit scheme which requires $\Delta t \leq O(\Delta S^2)$.

**Logarithmic change of variable** — $x = \ln S$ maps the geometric grid of stock prices onto a uniform arithmetic grid, making the diffusion coefficient constant and the discretisation simpler.

**Payoff polymorphism** — `Call` and `Put` inherit from an abstract `Payoff` class. Both the terminal condition and the boundary conditions are delegated to the payoff object, keeping the solver generic.
