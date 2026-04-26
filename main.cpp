#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <numeric> 

#include "Payoff.hpp"
#include "EDP.hpp"
#include "Res_full.hpp"
#include "Res_red.hpp"
#include "SDL.hpp"

// Fonction pour calculer l'erreur absolue entre deux vecteurs
std::vector<double> compute_error(const std::vector<double>& v1, const std::vector<double>& v2) {
    std::vector<double> err(v1.size());
    for(size_t i=0; i<v1.size(); ++i) {
        err[i] = std::abs(v1[i] - v2[i]);
    }
    return err;
}

// Fonction pour afficher les statistiques (Min, Max, Ecart-type)
void print_stats(const std::string& label, const std::vector<double>& data) {
    if (data.empty()) return;

    // Calcul Min et Max
    double min_val = *std::min_element(data.begin(), data.end());
    double max_val = *std::max_element(data.begin(), data.end());

    // Calcul Moyenne
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    double mean = sum / data.size();

    // Calcul Ecart-type
    double sq_sum = 0.0;
    for (double x : data) {
        sq_sum += (x - mean) * (x - mean);
    }
    double std_dev = std::sqrt(sq_sum / data.size());

    std::cout << "--- Statistiques " << label << " ---\n";
    std::cout << "  Min        : " << min_val << "\n";
    std::cout << "  Max        : " << max_val << "\n";
    std::cout << "  Moyenne    : " << mean << "\n";
    std::cout << "  Ecart-type : " << std_dev << "\n";
}

int main(int argc, char* argv[]) {
    (void)argc;
    (void)argv;
    
    const double T = 1.0, r = 0.1, sigma = 0.1, K = 100.0, L = 300.0;
    const int M = 1000, N = 1000;

    // Boucle sur les deux types d'options : Put et Call
    std::vector<std::unique_ptr<Payoff>> payoffs;
    payoffs.push_back(std::make_unique<Put>(K));
    payoffs.push_back(std::make_unique<Call>(K));

    for (auto& pay : payoffs) {
        std::string type = dynamic_cast<Put*>(pay.get()) ? "Put" : "Call";

        // Résolution Complète
        std::cout << "Calcul EDP complete..." << std::endl;
        EDP edp_full(T, r, sigma, L, M, N, *pay);
        ResolutionFull solver_full(edp_full);
        solver_full.solve();
        const auto& Cfull = edp_full.solution();
        const auto& S_grid = edp_full.S();

        //Résolution Réduite
        std::cout << "Calcul EDP reduite..." << std::endl;
        EDP edp_red(T, r, sigma, L, M, N, *pay);
        ResolutionReduced solver_red(edp_red);
        solver_red.solve();
        const auto& Cred = edp_red.solution();

        //Affichage console au strike K
        int idxK = static_cast<int>(K / (L / N));
        std::cout << "\nPrix au strike K (" << K << ") :\n";
        std::cout << "  C_full(0,K) = " << Cfull[idxK] << "\n";
        std::cout << "  C_red(0,K)  = " << Cred[idxK] << "\n\n";

        //Calcul de l'erreur absolue
        std::vector<double> diff = compute_error(Cfull, Cred);

        //Affichage des Statistiques
        print_stats("Methode Complete", Cfull);
        print_stats("Methode Reduite ", Cred);
        print_stats("Erreur Absolue " + type, diff);

        //Affichage SDL
        Sdl sdl;

        // Fenêtre 1 : courbes superposées
        sdl.w1 = std::make_unique<Window>("C(0,S) - " + type, 800, 600);
        sdl.w1->add_curve(S_grid, Cfull, 0, 0, 255);   // Bleu : complète
        sdl.w1->add_curve(S_grid, Cred, 255, 0, 0);    // Rouge : réduite

        // Fenêtre 2 : erreur absolue
        sdl.w2 = std::make_unique<Window>("Erreur absolue - " + type, 800, 600);
        sdl.w2->add_curve(S_grid, diff, 0, 255, 0);    // Vert

        // Affichage des fênetres
        sdl.show();
        sdl.wait_for_close();
    }

    return 0;
}