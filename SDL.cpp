#include "SDL.hpp"
#include <stdexcept>
#include <iostream>

/// @brief Vérifie si toutes les fenêtres sont fermées
static bool all_closed(const Sdl& s) {
    return !s.w1 && !s.w2;
}

/**
 * @brief Constructeur : initialise SDL.
 * @throws std::runtime_error si SDL_Init échoue
 */
Sdl::Sdl() {
    if (SDL_Init(SDL_INIT_VIDEO) != 0)
        throw std::runtime_error(SDL_GetError());
}

/**
 * @brief Destructeur : ferme toutes les fenêtres et quitte SDL.
 */
Sdl::~Sdl() {
    try { quit(); } catch (...) {}
}

/**
 * @brief Affiche le contenu de toutes les fenêtres actives.
 */
void Sdl::show() {
    if (w1) w1->render();
    if (w2) w2->render();
}

/**
 * @brief Ferme une fenêtre spécifique.
 * @param w Pointeur vers la fenêtre à fermer
 */
void Sdl::close(Window* w) {
    if (w1 && w1.get() == w) w1.reset();
    if (w2 && w2.get() == w) w2.reset();
}

/**
 * @brief Quitte SDL proprement.
 * 
 * Détruit toutes les fenêtres et appelle SDL_Quit().
 */
void Sdl::quit() {
    w1.reset(); 
    w2.reset();
    SDL_Quit();
    std::cout << "Fermeture de l'affichage..." << std::endl;
}

/**
 * @brief Boucle principale d'attente pour la fermeture des fenêtres.
 * 
 * Traite les événements SDL et met à jour l'affichage.
 */
void Sdl::wait_for_close() {
    SDL_Event e;
    while (!all_closed(*this)) {
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) {
                quit();
                return;
            }
            if (e.type == SDL_WINDOWEVENT &&
                e.window.event == SDL_WINDOWEVENT_CLOSE) {
                const Uint32 id = e.window.windowID;
                if (w1 && w1->id() == id) w1.reset();
                if (w2 && w2->id() == id) w2.reset();
            }
        }
        show();
        SDL_Delay(16); 
    }
}
