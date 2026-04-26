#ifndef THOMAS_HPP
#define THOMAS_HPP

#include <vector>

/**
 * @class ThomasAlgo
 * @brief Algorithme de Thomas pour résoudre un système tridiagonal.
 * 
 * Cette classe fournit une méthode statique pour résoudre des systèmes
 * linéaires tridiagonaux de la forme Ax = d, où A est tridiagonal.
 */
class ThomasAlgo {
public:
    /**
     * @brief Résout un système tridiagonal Ax = d avec l'algorithme de Thomas
     * @param a Diagonale inférieure
     * @param b Diagonale principale
     * @param c Diagonale supérieure
     * @param d Membre de droite 
     * @return Le vecteur solution x
     */
    static std::vector<double> solve(
        const std::vector<double>& a, 
        const std::vector<double>& b, 
        const std::vector<double>& c, 
        const std::vector<double>& d
    );
};

#endif