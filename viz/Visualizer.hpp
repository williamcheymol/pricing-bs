#pragma once
#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <vector>
#include <string>

// =============================================================================
//  Visualizer.hpp  —  Real-time SDL2 window for BS-Pricer
// =============================================================================
//
//  Six display modes (keys 1–6):
//    1. European Call      — FD curve + BS analytical overlay
//    2. European Put       — FD curve + BS analytical overlay
//    3. Digital Call       — FD + Digital Put overlay (parity check)
//    4. American Put       — FD + European Put overlay (early-exercise premium)
//    5. Up-and-Out Call    — FD + Vanilla Call overlay, barrier line
//    6. Down-and-Out Put   — FD + Vanilla Put overlay, barrier line
//
//  Interactive controls:
//    ↑ / ↓     sigma  ± 0.01   (range  0.01 – 0.80)
//    ← / →     r      ± 0.005  (range  0.001 – 0.20)
//    W / S     T      ± 0.1    (range  0.1  – 5.0)
//    1 – 6     switch display mode
//    ESC / Q   quit
// =============================================================================

// CurveData: a set of (S, V) pairs for one plotted curve.
struct CurveData {
    std::vector<double> S;
    std::vector<double> V;
};

class Visualizer {
public:
    Visualizer();
    ~Visualizer();

    bool init();   // open window, renderer, load font; returns false on failure
    void run();    // enter event loop (blocks until window is closed)

private:
    // ----- SDL resources -----
    SDL_Window*   win_    = nullptr;
    SDL_Renderer* ren_    = nullptr;
    TTF_Font*     font_   = nullptr;   // ~16 pt, titles + values
    TTF_Font*     fontsm_ = nullptr;   // ~13 pt, axis labels + help text

    // ----- interactive model parameters -----
    double sigma_ = 0.20;
    double r_     = 0.05;
    double T_     = 1.00;
    double K_     = 100.0;   // A/D to adjust ±5
    int    mode_  = 0;   // 0–5

    // Solver toggle. R to switch between CN and ReducedCN.
    bool use_rcn_ = false;

    // Payoff smoothing toggle (Digital modes only). G to toggle.
    bool smooth_ = true;

    // Curve visibility toggles. M = main, O = overlay.
    bool show_main_    = true;
    bool show_overlay_ = true;

    // ----- solver grid resolution (speed vs accuracy trade-off) -----
    static constexpr int    NV    = 300;
    static constexpr int    MV    = 300;
    static constexpr double SMAX  = 250.0;

    // ----- window / plot layout (pixels) -----
    static constexpr int WIN_W = 1050;
    static constexpr int WIN_H = 680;
    static constexpr int PX0   = 80;    // plot left edge
    static constexpr int PY0   = 52;    // plot top edge
    static constexpr int PW    = 650;   // plot width
    static constexpr int PH    = 500;   // plot height

    // ----- state -----
    bool running_ = false;
    bool dirty_   = true;   // true → re-solve next frame (expensive)

    int mouse_x_ = -1;      // current mouse position (-1 = outside window)
    int mouse_y_ = -1;

    // ----- scene (built once per dirty frame, re-used until next dirty) -----
    struct Scene {
        CurveData  main, overlay, intrinsic;
        double     S_lo = 0.0, S_hi = 200.0;
        double     V_lo = 0.0, V_hi = 50.0;
        double     barrier = -1.0;   // <0 = no barrier line
        std::string title;
        std::string main_label;
        std::string overlay_label;
        bool has_overlay   = false;
        bool has_intrinsic = false;
    };

    // ----- internal methods -----
    Scene  scene_;   // cached — rebuilt only when dirty_, redrawn every frame

    void   handle_events();
    void   render(const Scene& sc);
    Scene  build_scene();
    void   draw_hover(const Scene& sc);

    // drawing primitives
    void fill_rect(int x, int y, int w, int h, SDL_Color c);
    void draw_line(int x1, int y1, int x2, int y2, SDL_Color c);
    void draw_curve(const CurveData& cd, const Scene& sc,
                    SDL_Color col, int thickness = 1, bool dashed = false);
    void draw_vline_dashed(double S, const Scene& sc, SDL_Color col);
    void draw_axes(const Scene& sc);
    void draw_side_panel(const Scene& sc);
    void draw_legend(const Scene& sc);
    void draw_help();
    void draw_text(const std::string& s, int x, int y, SDL_Color col,
                   bool large = false);

    // world → pixel coordinate transforms
    int sx(double S, const Scene& sc) const;   // S → screen x
    int vx(double V, const Scene& sc) const;   // V → screen y (inverted)

    // font loader: tries common system paths; returns nullptr if all fail
    static TTF_Font* try_load_font(int size);
};
