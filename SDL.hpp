#ifndef SDL_HPP
#define SDL_HPP

#include <memory>
#include "SDL2/SDL.h"
#include "window.hpp"
#include <iostream>


/**
 * @class Sdl
 * @brief Gestionnaire principal de l'environnement graphique SDL.
 * * Cette classe s'occupe de l'initialisation de la bibliothèque SDL2 et de la gestion de la mémoire des fenêtres.
 */
class Sdl {
public:
    /** 
     * @brief Première fenêtre : Comparaison graphique.
     * @details Affiche la superposition des méthodes complète et réduite .
    */
    std::unique_ptr<Window> w1;

    /** 
     * @brief Deuxième fenêtre : Analyse de l'erreur.
     * @details Affiche la différence absolue entre les deux méthodes .
    */
    std::unique_ptr<Window> w2;

    /**
     * @brief Constructeur.
     * * Initialise le sous-système vidéo de la SDL .
     * Lance une exception  si l'initialisation échoue.
     */
    Sdl();

    /**
     * @brief Destructeur.
     * * Ferme automatiquement toutes les fenêtres et quitte la SDL.
     */
    ~Sdl();
    
    /**
     * @brief Affiche le contenu de toutes les fenêtres actives.
    */
    void show();
    
    /**
     * @brief Boucle d'attente principale (Event Loop).
     * * Bloque le programme et attend que l'utilisateur ferme toutes les fenêtres.
    */    
    void wait_for_close();

    /**
     * @brief Ferme une fenêtre spécifique.
     * * Recherche quel unique_ptr possède le pointeur w et le réinitialise.
     * @param w Pointeur b vers la fenêtre à fermer.
    */
    void close(Window* w);     

    /**
     * @brief Quitte proprement la SDL.
     * * Détruit toutes les fenêtres et appelle SDL_Quit().
    */
    void quit();               
};

#endif
