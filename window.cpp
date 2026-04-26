#include "window.hpp"
#include <iostream>
#include <limits>

/**
 * @brief Constructeur de Window.
 * 
 * Crée une fenêtre SDL et un renderer associé.
 * Affiche un message d'erreur si la création échoue.
 */
Window::Window(const std::string& title, int w, int h)
    : width_(w), height_(h)
{
    window_ = SDL_CreateWindow(title.c_str(),
                               SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                               width_, height_,
                               SDL_WINDOW_SHOWN);
    if (!window_) {
        std::cerr << "Erreur Window creation: " << SDL_GetError() << std::endl;
    }

    renderer_ = SDL_CreateRenderer(window_, -1, SDL_RENDERER_ACCELERATED);
    if (!renderer_) {
        std::cerr << "Erreur Renderer creation: " << SDL_GetError() << std::endl;
    }
}

/**
 * @brief Destructeur.
 * 
 * Libère le renderer et la fenêtre.
 */
Window::~Window() {
    if (renderer_) SDL_DestroyRenderer(renderer_);
    if (window_) SDL_DestroyWindow(window_);
}

/**
 * @brief Récupère l'identifiant SDL de la fenêtre.
 * @return L'ID SDL de la fenêtre.
 */
Uint32 Window::id() const {
    return SDL_GetWindowID(window_);
}

/**
 * @brief Ajoute une courbe à tracer.
 */
void Window::add_curve(const std::vector<double>& x, const std::vector<double>& y, Uint8 r, Uint8 g, Uint8 b) {
    if (x.size() != y.size()) {
        std::cerr << "Erreur taille vecteurs diff (Window::add_curve)" << std::endl;
        return;
    }
    curves_.push_back({x, y, r, g, b});
}

/**
 * @brief Calcule les bornes globales des courbes.
 */
void Window::get_min_max(double& x_min, double& x_max, double& y_min, double& y_max) const {
    x_min = std::numeric_limits<double>::max();
    x_max = std::numeric_limits<double>::lowest();
    y_min = std::numeric_limits<double>::max();
    y_max = std::numeric_limits<double>::lowest();

    if (curves_.empty()) {
        x_min = 0; x_max = 1; y_min = 0; y_max = 1;
        return;
    }

    for (const auto& c : curves_) {
        for (double val : c.x) {
            if (val < x_min) x_min = val;
            if (val > x_max) x_max = val;
        }
        for (double val : c.y) {
            if (val < y_min) y_min = val;
            if (val > y_max) y_max = val;
        }
    }
}

/**
 * @brief Convertit des coordonnées mathématiques en pixels.
 */
void Window::map_to_screen(double x, double y,
                           double x_min, double x_max, double y_min, double y_max,
                           int& px, int& py) const 
{
    int mx = width_ * 0.1;
    int my = height_ * 0.1;
    int draw_w = width_ - 2 * mx;
    int draw_h = height_ - 2 * my;

    double rx = (x - x_min) / (x_max - x_min);
    double ry = (y - y_min) / (y_max - y_min);

    px = mx + static_cast<int>(rx * draw_w);
    py = (height_ - my) - static_cast<int>(ry * draw_h);
}

/**
 * @brief Affiche toutes les courbes sur la fenêtre.
 */
void Window::render() {
    if (!renderer_) return;

    // Efface l'écran en blanc
    SDL_SetRenderDrawColor(renderer_, 255, 255, 255, 255);
    SDL_RenderClear(renderer_);

    // Bornes des courbes
    double x_min, x_max, y_min, y_max;
    get_min_max(x_min, x_max, y_min, y_max);

    // Dessine axes
    SDL_SetRenderDrawColor(renderer_, 0, 0, 0, 255);
    int mx = width_ * 0.1;
    int my = height_ * 0.1;
    SDL_RenderDrawLine(renderer_, mx, my, mx, height_ - my);
    SDL_RenderDrawLine(renderer_, mx, height_ - my, width_ - mx, height_ - my);

    // Dessine les courbes
    for (const auto& c : curves_) {
        SDL_SetRenderDrawColor(renderer_, c.r, c.g, c.b, 255);
        if (c.x.size() < 2) continue;

        int px1, py1, px2, py2;
        map_to_screen(c.x[0], c.y[0], x_min, x_max, y_min, y_max, px1, py1);

        for (size_t i = 1; i < c.x.size(); ++i) {
            map_to_screen(c.x[i], c.y[i], x_min, x_max, y_min, y_max, px2, py2);
            SDL_RenderDrawLine(renderer_, px1, py1, px2, py2);
            px1 = px2;
            py1 = py2;
        }
    }

    SDL_RenderPresent(renderer_);
}
