#include "Res_full.hpp"

/**
 * @brief Constructeur de ResolutionFull.
 * @param edp Référence vers l'EDP à résoudre
 */
ResolutionFull::ResolutionFull(EDP& edp) : edp_(edp) {}

//Fonctions auxiliaires

double ResolutionFull::A(int i) const {
    const double S = i * edp_.dS();
    const double sig = edp_.sigma();
    const double r = edp_.r();
    const double dt = edp_.dt();
    const double dS = edp_.dS();
    return 0.25 * dt * ((sig*sig*S*S)/(dS*dS) - (r*S)/dS);
}

double ResolutionFull::B(int i) const {
    const double S = i * edp_.dS();
    const double sig = edp_.sigma();
    const double r = edp_.r();
    const double dt = edp_.dt();
    const double dS = edp_.dS();
    return -0.5 * dt * ((sig*sig*S*S)/(dS*dS) + r);
}

double ResolutionFull::C(int i) const {
    const double S = i * edp_.dS();
    const double sig = edp_.sigma();
    const double r = edp_.r();
    const double dt = edp_.dt();
    const double dS = edp_.dS();
    return 0.25 * dt * ((sig*sig*S*S)/(dS*dS) + (r*S)/dS);
}

/**
 * @brief Construit le vecteur du membre de droite du système tridiagonal
 *        pour un pas de temps donné, incluant les contributions des bords.
 * @param N Nombre de points spatiaux
 * @param C0_np1 Bord gauche au temps n+1
 * @param CN_np1 Bord droit au temps n+1
 * @param C0_n Bord gauche au temps n
 * @param CN_n Bord droit au temps n
 * @return Vecteur du membre de droite
 */
std::vector<double> ResolutionFull::buildRHS(int N, double C0_np1, double CN_np1, double C0_n, double CN_n) const {
    std::vector<double> d(N-1);
    const auto& C_next = edp_.C(); // Solution au temps n+1

    for (int i = 1; i <= N-1; ++i) {
        const int k = i-1;
        double Cim1 = (i == 1) ? C0_np1 : C_next[i-1];
        double Ci   = C_next[i];
        double Cip1 = (i == N-1) ? CN_np1 : C_next[i+1];
        d[k] = A(i)*Cim1 + (1.0 + B(i))*Ci + C(i)*Cip1;
    }

    // Contribution des bords au temps n
    d[0]   += A(1) * C0_n;
    d[N-2] += C(N-1) * CN_n;

    return d;
}

/**
 * @brief Construit les vecteurs tridiagonaux a,b,c pour un pas de temps donné.
 * @param N Nombre de points spatiaux
 * @param a Vecteur des coefficients inférieurs
 * @param b Vecteur des coefficients diagonaux
 * @param c Vecteur des coefficients supérieurs
 */
void ResolutionFull::buildTridiagonal(int N, std::vector<double>& a, std::vector<double>& b, std::vector<double>& c) const {
    a.resize(N-1); b.resize(N-1); c.resize(N-1);
    for (int i = 1; i <= N-1; ++i) {
        int k = i-1;
        a[k] = -A(i);
        b[k] = 1.0 - B(i);
        c[k] = -C(i);
    }
}

/**
 * @brief Résolution complète de l'EDP sur toute la grille temporelle.
 * 
 * Utilise la méthode de Crank--Nicolson et l'algorithme de Thomas
 * pour résoudre les systèmes tridiagonaux à chaque pas de temps.
 * La solution finale est stockée dans l'objet EDP.
 */
void ResolutionFull::solve() {
    const int N = edp_.N();
    const int M = edp_.M();

    std::vector<double> a, b, c, d;

    for (int n = M-1; n >= 0; --n) {
        const double t_np1 = (n+1) * edp_.dt();
        const double t_n   = n * edp_.dt();

        const double C0_np1 = edp_.payoff().price_zero(t_np1, edp_.T(), edp_.r());
        const double CN_np1 = edp_.payoff().price_limit(t_np1, edp_.T(), edp_.r());
        const double C0_n   = edp_.payoff().price_zero(t_n, edp_.T(), edp_.r());
        const double CN_n   = edp_.payoff().price_limit(t_n, edp_.T(), edp_.r());

        buildTridiagonal(N, a, b, c);
        d = buildRHS(N, C0_np1, CN_np1, C0_n, CN_n);

        // Résolution du système tridiagonal
        std::vector<double> U_n = ThomasAlgo::solve(a, b, c, d);

        // Mise à jour de la solution
        edp_.update_interior(U_n, t_n, edp_.r());
    }
}
