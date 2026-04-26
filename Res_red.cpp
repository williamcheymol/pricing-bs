#include "Res_red.hpp"

/**
 * @brief Constructeur du solveur réduit.
 * @param edp Référence vers l'objet EDP à résoudre
 */
ResolutionReduced::ResolutionReduced(EDP& edp)
    : edp_(edp) {}

/**
 * @brief Résout l'équation de Black--Scholes via l'EDP réduite.
 * 
 * Étapes :
 * 1. Calcul des coefficients (alpha, beta, lambda...)
 * 2. Initialisation (grille logarithmique et condition initiale)
 * 3. Résolution temporelle avec Thomas
 * 4. Transformation inverse et interpolation sur la grille originale
 */
void ResolutionReduced::solve() {
    // 1. Calcul des coefficients
    HeatParams p = compute_params();

    // 2. Initialisation (grille x et condition initiale u)
    std::vector<double> u(edp_.N() + 1);
    std::vector<double> x_grid(edp_.N() + 1);
    initialize_grid(p, u, x_grid);

    // 3. Résolution temporelle (schéma implicite + Thomas)
    perform_time_stepping(p, u);

    // 4. Transformation inverse vers la grille originale et interpolation
    transform_and_interpolate(p, u, x_grid);
}

/**
 * @brief Calcule tous les paramètres nécessaires au changement de variable logarithmique.
 * @return Structure HeatParams remplie
 */
ResolutionReduced::HeatParams ResolutionReduced::compute_params() const {
    HeatParams p;
    p.K = 100.0; // Strike pivot
    double sigma2 = edp_.sigma() * edp_.sigma();

    // Paramètres temporels
    p.tau_final = 0.5 * sigma2 * edp_.T();
    double dtau = p.tau_final / edp_.M();

    // Paramètres spatiaux (log)
    double S_min_log = 0.1;
    p.x_min = std::log(S_min_log / p.K);
    double x_max = std::log(edp_.L() / p.K);
    p.dx = (x_max - p.x_min) / edp_.N();

    // Paramètre pour le schéma implicite
    p.lambda = dtau / (p.dx * p.dx);

    // Coefficients alpha et beta
    double k_param = 2.0 * edp_.r() / sigma2;
    p.alpha = -(k_param - 1.0) / 2.0;
    p.beta = -std::pow(k_param + 1.0, 2) / 4.0;

    return p;
}

/**
 * @brief Initialise la grille logarithmique et la condition initiale transformée u(0,x).
 * @param[in] p Paramètres calculés
 * @param[out] u Vecteur solution à initialiser
 * @param[out] x_grid Grille spatiale logarithmique
 */
void ResolutionReduced::initialize_grid(const HeatParams& p, 
                                        std::vector<double>& u, 
                                        std::vector<double>& x_grid) const {
    for (int i = 0; i <= edp_.N(); ++i) {
        x_grid[i] = p.x_min + i * p.dx;
        double S_val = p.K * std::exp(x_grid[i]);
        double val_payoff = edp_.payoff()(S_val);
        u[i] = (val_payoff / p.K) * std::exp(-p.alpha * x_grid[i]);
    }
}

/**
 * @brief Boucle temporelle avec le schéma implicite. Utilise l'algorithme de Thomas.
 * @param[in] p Paramètres de simulation
 * @param[in,out] u Solution mise à jour à chaque pas de temps
 */
void ResolutionReduced::perform_time_stepping(const HeatParams& p, std::vector<double>& u) const {
    int N = edp_.N();
    int M = edp_.M();

    std::vector<double> a(N - 1), b(N - 1), c(N - 1), d(N - 1);

    for (int m = 1; m <= M; ++m) {
        for (int i = 0; i < N - 1; ++i) {
            a[i] = -p.lambda;
            b[i] = 1.0 + 2.0 * p.lambda;
            c[i] = -p.lambda;
            d[i] = u[i + 1];
        }

        // Conditions aux bords
        d[0] -= a[0] * u[0];
        d[N - 2] -= c[N - 2] * u[N];

        // Résolution du système tridiagonal
        std::vector<double> sol = ThomasAlgo::solve(a, b, c, d);

        for (int i = 0; i < N - 1; ++i) {
            u[i + 1] = sol[i];
        }
    }
}

/**
 * @brief Transforme la solution logarithmique u(x) vers la grille originale et interpoler.
 * @param[in] p Paramètres du changement de variable
 * @param[in] u Solution finale sur la grille logarithmique
 * @param[in] x_grid Grille logarithmique correspondante
 */
void ResolutionReduced::transform_and_interpolate(const HeatParams& p, 
                                                  const std::vector<double>& u, 
                                                  const std::vector<double>& x_grid) {
    int N = edp_.N();
    double time_factor = std::exp(p.beta * p.tau_final);
    std::vector<double> C_log(N + 1);

    for (int i = 0; i <= N; ++i) {
        C_log[i] = p.K * u[i] * std::exp(p.alpha * x_grid[i]) * time_factor;
    }

    // Interpolation sur la grille originale
    for (int i = 0; i <= N; ++i) {
        double target_S = edp_.S()[i];

        if (target_S <= std::exp(p.x_min) * p.K) {
            edp_.C()[i] = edp_.payoff().price_zero(0, edp_.T(), edp_.r());
        } else if (target_S >= edp_.L()) {
            edp_.C()[i] = C_log[N];
        } else {
            double target_x = std::log(target_S / p.K);
            double pos = (target_x - p.x_min) / p.dx;
            int idx = static_cast<int>(pos);
            double w = pos - idx;
            edp_.C()[i] = (1.0 - w) * C_log[idx] + w * C_log[idx + 1];
        }
    }
}
