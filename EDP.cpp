#include "EDP.hpp"
#include "Res_full.hpp"
#include <stdexcept>

/**
 * @brief Constructeur de la classe EDP.
 * 
 * Initialise la discrétisation en temps et en espace et calcule la grille
 * des prix sous-jacents S_[i].
 * Initialise également l'état terminal via initTerminal().
 */
EDP::EDP(double T, double r, double sigma, double L,
         int M, int N, const Payoff& payoff)
    : T_(T), r_(r), sigma_(sigma), L_(L),
      M_(M), N_(N), payoff_(payoff)
{
    dt_ = T_ / static_cast<double>(M_);
    dS_ = L_ / static_cast<double>(N_);

    S_.resize(N_ + 1);
    C_.resize(N_ + 1);

    // Remplissage de la grille spatiale
    for (int i = 0; i <= N_; ++i) {
        S_[i] = i * dS_;
    }

    initTerminal(); // Initialise C(T,S)
}

/**
 * @brief Initialise la condition terminale C(T,S) sur la grille.
 */
void EDP::initTerminal() {
    for (int i = 0; i <= N_; ++i) {
        C_[i] = payoff_(S_[i]);
    }
}

/**
 * @brief Impose les conditions aux bords pour un temps donné.
 * @param t Temps actuel
 * @param r Taux sans risque
 */
void EDP::applyBoundary(double t, double r) {
    C_[0]   = payoff_.price_zero(t, T_, r);   // Bord S=0
    C_[N_]  = payoff_.price_limit(t, T_, r);  // Bord S=L
}

/**
 * @brief Met à jour uniquement les points intérieurs de la grille
 * et applique les conditions aux bords.
 * @param interior Vecteur des nouvelles valeurs intérieures
 * @param t Temps actuel
 * @param r Taux sans risque
 * @throws std::invalid_argument si la taille de interior est incorrecte
 */
void EDP::update_interior(const std::vector<double>& interior, double t, double r) {
    if (static_cast<int>(interior.size()) != N_ - 1) {
        throw std::invalid_argument("update_interior(): taille incorrecte");
    }

    for (int i = 1; i <= N_ - 1; ++i) {
        C_[i] = interior[i - 1];
    }

    applyBoundary(t, r);
}

/**
 * @brief Met à jour toute la grille avec next et applique les conditions aux bords.
 * @param next Vecteur de la nouvelle solution
 * @param t Temps actuel
 * @param r Taux sans risque
 * @throws std::invalid_argument si la taille de next est incorrecte
 */
void EDP::update(const std::vector<double>& next, double t, double r) {
    if (static_cast<int>(next.size()) != N_ + 1) {
        throw std::invalid_argument("update(): taille incorrecte");
    }

    C_ = next;
    applyBoundary(t, r);
}

/**
 * @brief Retourne la solution actuelle C(t,S)
 * @return Référence constante vers le vecteur de solution
 */
const std::vector<double>& EDP::solution() const {
    return C_;
}