#ifndef RES_FULL_HPP
#define RES_FULL_HPP

#include "EDP.hpp"
#include "Thomas.hpp"
#include <vector>

/**
 * @class ResolutionFull
 * @brief Classe responsable de la résolution de la forme complète de l'EDP.
 * 
 * Cette classe utilise l'objet EDP fourni pour appliquer le schéma de Crank--Nicolson et calculer la solution de l'équation de Black--Scholes sur toute la grille.
 * Elle utilise un algorithme tridiagonal (Thomas) pour résoudre les systèmes linéaires à chaque pas de temps.
 */
class ResolutionFull {
public:
    /**
     * @brief Constructeur de ResolutionFull.
     * @param edp Référence vers l'EDP à résoudre.
     */
    ResolutionFull(EDP& edp);

    /**
     * @brief Résout l'EDP sur toute la grille temporelle et spatiale.
     * 
     * La solution finale est stockée dans l'objet EDP fourni.
     */
    void solve();

private:
    /// Référence à l'EDP à résoudre
    EDP& edp_;

    //Fonctions auxiliaires pour le calcul des coefficients
    double A(int i) const;
    double B(int i) const;
    double C(int i) const;

    /**
     * @brief Construit les vecteurs de la matrice tridiagonale a,b,c.
     */
    void buildTridiagonal(int N, std::vector<double>& a, std::vector<double>& b, std::vector<double>& c) const;

    /**
     * @brief Construit le vecteur du membre de droite d du système tridiagonal.
     */
    std::vector<double> buildRHS(int N, double C0_np1, double CN_np1, double C0_n, double CN_n) const;
};

#endif
