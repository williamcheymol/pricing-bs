#include "Thomas.hpp"

/**
 * @brief Résout un système tridiagonal Ax = d avec l'algorithme de Thomas.
 * 
 * Implémente la méthode de factorisation LU simplifiée pour les matrices tridiagonales.
 * @param a Diagonale inférieure
 * @param b Diagonale principale
 * @param c Diagonale supérieure
 * @param d Membre de droite
 * @return Vecteur solution x
 */
std::vector<double> ThomasAlgo::solve(
    const std::vector<double>& a, 
    const std::vector<double>& b, 
    const std::vector<double>& c, 
    const std::vector<double>& d
) {
    std::vector<double> c_prime = c;
    std::vector<double> d_prime = d;

    int n = d_prime.size();
    int m = n - 1; 

    // Première ligne
    c_prime[0] /= b[0];
    d_prime[0] /= b[0];

    // Forward sweep
    for (int i = 1; i < m; i++) {
        double temp = b[i] - a[i] * c_prime[i - 1];
        c_prime[i] /= temp;
        d_prime[i] = (d_prime[i] - a[i] * d_prime[i - 1]) / temp;
    }

    // Dernière ligne
    d_prime[m] = (d_prime[m] - a[m] * d_prime[m - 1]) / (b[m] - a[m] * c_prime[m - 1]);

    // Back substitution
    for (int i = m - 1; i >= 0; i--) {
        d_prime[i] -= c_prime[i] * d_prime[i + 1];
    }

    return d_prime;
}
