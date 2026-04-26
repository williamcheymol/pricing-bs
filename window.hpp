#ifndef WINDOW_HPP
#define WINDOW_HPP

#include <SDL2/SDL.h>
#include <vector>
#include <string>
#include <algorithm>

/**
 * @class Window
 * @brief Classe encapsulant une fenêtre SDL2 pour le tracé de courbes 2D.
 * * Cette classe gère la création et la destruction de la fenêtre et affiche plusieurs courbes superposées.
 */
class Window {
public:
    /**
     * @brief Constructeur.
     * * Initialise la fenêtre SDL et le moteur de rendu (Renderer).
     * @param title Le titre de la fenêtre.
     * @param width Largeur de la fenêtre en pixels.
     * @param height Hauteur de la fenêtre en pixels.
    */
    Window(const std::string& title, int width, int height);
    
    /**
     * @brief Destructeur.
     * * Libère la mémoire SDL .
    */
    ~Window();

    /**
     * @brief Dessine le contenu de la fenêtre.
    */ 
    void render();

    /**
     * @brief Récupère l'identifiant unique SDL de la fenêtre.
     * * Nécessaire pour savoir quelle fenêtre l'utilisateur veut fermer.
     * @return L'ID  de la fenêtre.
    */    
    Uint32 id() const;

    /**
     * @brief Ajoute une courbe à la liste des éléments à tracer.
     * @param x Vecteur des abscisses.
     * @param y Vecteur des ordonnées.
     * @param r Composante Rouge de la couleur (0-255).
     * @param g Composante Verte de la couleur (0-255).
     * @param b Composante Bleue de la couleur (0-255).
    */
    void add_curve(const std::vector<double>& x, const std::vector<double>& y, 
                   Uint8 r, Uint8 g, Uint8 b);

private:
    SDL_Window* window_;
    SDL_Renderer* renderer_;
    int width_;
    int height_;

    /**
     * @struct Curve
     * @brief Structure interne pour stocker les données et la couleur d'une courbe.
    */
    struct Curve {
        std::vector<double> x;
        std::vector<double> y;
        Uint8 r, g, b;
    };
    std::vector<Curve> curves_;

    /**
     * @brief Calcule les bornes  globales de toutes les courbes.
     * * Permet d'adapter l'échelle du graphique.
     * @param[out] x_min Minimum des abscisses trouvé.
     * @param[out] x_max Maximum des abscisses trouvé.
     * @param[out] y_min Minimum des ordonnées trouvé.
     * @param[out] y_max Maximum des ordonnées trouvé.
    */
    void get_min_max(double& x_min, double& x_max, double& y_min, double& y_max) const;

    /**
     * @brief Convertit des coordonnées mathématiques (x,y) en coordonnées pixels (px,py).
     * @param x Abscisse mathématique à convertir.
     * @param y Ordonnée mathématique à convertir.
     * @param x_min Borne min des abscisses.
     * @param x_max Borne max des abscisses.
     * @param y_min Borne min des ordonnées.
     * @param y_max Borne max des ordonnées.
     * @param[out] px Coordonnée X en pixels résultante.
     * @param[out] py Coordonnée Y en pixels résultante.
     */
    void map_to_screen(double x, double y, 
                       double x_min, double x_max, double y_min, double y_max,
                       int& px, int& py) const;
};

#endif

