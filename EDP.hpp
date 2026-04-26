#ifndef EDP_HPP
#define EDP_HPP

#include <vector>
#include "Payoff.hpp"

/**
 * @class EDP
 * @brief Classe représentant une équation aux dérivées partielles (EDP) pour le pricing d'options.
 * 
 * Cette classe gère la discrétisation spatio-temporelle de l'équation de
 * Black-Scholes et stocke la solution intermédiaire sur la grille.
 * Elle dépend d'un objet Payoff pour définir les conditions terminales et aux bords.
 */
class EDP {
protected:
    /// Temps terminal de l'option
    double T_;
    /// Taux d'intérêt sans risque
    double r_;
    /// Volatilité de l'actif sous-jacent
    double sigma_;
    /// Borne maximale du prix de l'actif 
    double L_;
    /// Nombre de pas de temps
    int M_;
    /// Nombre de pas en espace
    int N_;

    /// Pas de temps
    double dt_;
    /// Pas spatial
    double dS_;

    /// Grille spatiale des prix de l'actif
    std::vector<double> S_;
    /// Grille de la solution C(t,S)
    std::vector<double> C_;

    /// Payoff associé à l'option
    const Payoff& payoff_;

public:
    /**
     * @brief Constructeur de la classe EDP.
     * @param T Temps terminal
     * @param r Taux d'intérêt sans risque
     * @param sigma Volatilité
     * @param L Borne maximale du prix du sous-jacent
     * @param M Nombre de pas de temps
     * @param N Nombre de pas en espace
     * @param payoff Référence au payoff associé
     */
    EDP(double T, double r, double sigma, double L, int M, int N, const Payoff& payoff);

    /// Destructeur virtuel par défaut
    virtual ~EDP() = default;

    /**
     * @brief Retourne la solution numérique C(0,S).
     * @return Référence constante vers le vecteur des valeurs de l'option.
     */
    const std::vector<double>& solution() const;

    // Getters scalaires

    /**
     * @brief Retourne la maturité de l'option.
     * @return Temps final T.
     */
    double T() const { return T_; }

    /**
     * @brief Retourne la borne maximale du sous-jacent.
     * @return Valeur maximale L.
     */
    double L() const { return L_; }

    /**
     * @brief Retourne le nombre de pas de temps.
     * @return Nombre de pas temporels M.
     */
    int M() const { return M_; }

    /**
     * @brief Retourne le nombre de pas d'espace.
     * @return Nombre de pas spatiaux N.
     */
    int N() const { return N_; }

    /**
     * @brief Retourne le pas de temps.
     * @return Valeur de dt.
     */
    double dt() const { return dt_; }

    /**
     * @brief Retourne le pas d'espace.
     * @return Valeur de dS.
     */
    double dS() const { return dS_; }

    /**
     * @brief Retourne le taux d'intérêt sans risque.
     * @return Taux r.
     */
    double r() const { return r_; }

    /**
     * @brief Retourne la volatilité du sous-jacent.
     * @return Volatilité sigma.
     */
    double sigma() const { return sigma_; }

    /**
     * @brief Retourne le payoff associé à l'EDP.
     * @return Référence constante vers le payoff.
     */
    const Payoff& payoff() const { return payoff_; }

    // Getters pour les grilles

    /**
     * @brief Retourne la grille spatiale du sous-jacent.
     * @return Référence constante vers le vecteur S.
     */
    const std::vector<double>& S() const { return S_; }

    /**
     * @brief Retourne la grille de solution C(0,S).
     * @return Référence constante vers le vecteur C.
     */
    const std::vector<double>& C() const { return C_; }

    /**
     * @brief Accès modifiable à la grille de solution.
     * @return Référence vers le vecteur C.
     */
    std::vector<double>& C() { return C_; }

    /**
     * @brief Initialise la condition terminale C(T,S) sur la grille
     */
    void initTerminal();

    /**
     * @brief Impose les conditions aux bords à un temps donné
     * @param t Temps actuel
     * @param r Taux sans risque
     */
    void applyBoundary(double t, double r);

    /**
     * @brief Met à jour uniquement les points intérieurs de la grille
     * @param interior Vecteur des nouvelles valeurs intérieures
     * @param t Temps actuel
     * @param r Taux sans risque
     */
    void update_interior(const std::vector<double>& interior, double t, double r);

    /**
     * @brief Mise à jour : remplace la grille courante par next et impose les bords
     * @param next Vecteur de la nouvelle solution
     * @param t Temps actuel
     * @param r Taux sans risque
     */
    void update(const std::vector<double>& next, double t, double r);
};

#endif
