#ifndef RES_RED_HPP
#define RES_RED_HPP

#include "EDP.hpp"
#include "Thomas.hpp"
#include <vector>
#include <cmath>

/**
 * @class ResolutionReduced
 * @brief Solveur de l'équation de Black--Scholes utilisant la méthode de l'EDP réduite.
 * 
 * Cette classe possède une référence vers un objet EDP et applique un changement
 * de variable logarithmique pour simplifier l'équation.
 * Elle utilise l'algorithme tridiagonal de Thomas pour résoudre la matrice à chaque pas de temps.
 */
class ResolutionReduced {
public:
    /**
     * @brief Constructeur du solveur réduit.
     * @param edp Référence vers l'objet EDP à résoudre.
     */
    ResolutionReduced(EDP& edp);

    /**
     * @brief Résout l'équation de Black-Scholes avec l'EDP réduite.
     */
    void solve();

private:
    /// Référence à l'EDP utilisé
    EDP& edp_;

    /**
     * @struct HeatParams
     * @brief Paramètres internes pour le changement de variable logarithmique.
     */
    struct HeatParams {
        double alpha, beta, lambda;
        double tau_final, dx, x_min;
        double K;
    };

    /**
     * @brief Calcule les paramètres nécessaires au changement de variable.
     * @return Une structure HeatParams remplie.
     */
    HeatParams compute_params() const;

    /**
     * @brief Initialise la grille spatiale et la condition initiale transformée.
     * @param[in] p Paramètres calculés.
     * @param[out] u Vecteur solution à initialiser.
     * @param[out] x_grid Grille spatiale logarithmique.
     */
    void initialize_grid(const HeatParams& p, 
                         std::vector<double>& u, 
                         std::vector<double>& x_grid) const;

    /**
     * @brief Effectue la boucle temporelle avec le schéma implicite.
     * Utilise l'algorithme de Thomas pour inverser la matrice tridiagonale à chaque étape.
     * @param[in] p Paramètres de simulation.
     * @param[in,out] u Vecteur solution mis à jour à chaque pas de temps.
     */
    void perform_time_stepping(const HeatParams& p, 
                               std::vector<double>& u) const;

    /**
     * @brief Effectue le changement de variable inverse et l'interpolation sur la grille d'origine.
     * @param[in] p Paramètres de transformation.
     * @param[in] u Solution finale de l'équation.
     * @param[in] x_grid Grille spatiale logarithmique.
     */
    void transform_and_interpolate(const HeatParams& p, 
                                   const std::vector<double>& u, 
                                   const std::vector<double>& x_grid);
};

#endif
