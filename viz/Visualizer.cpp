#include "Visualizer.hpp"
#include "Option.hpp"
#include "Grid.hpp"
#include "CrankNicolson.hpp"
#include "ReducedCN.hpp"
#include "Solver.hpp"
#include "BSFormula.hpp"

#include <cmath>
#include <algorithm>
#include <limits>
#include <sstream>
#include <iomanip>

// =============================================================================
//  Colour palette  (R, G, B, A)
// =============================================================================
static constexpr SDL_Color C_BG       = {13,  17,  23, 255};   // dark navy
static constexpr SDL_Color C_PLOTBG   = {18,  26,  44, 255};   // slightly lighter
static constexpr SDL_Color C_GRIDLN   = {32,  48,  72, 255};   // faint blue grid
static constexpr SDL_Color C_PANEL    = {20,  30,  52, 255};   // side-panel bg
static constexpr SDL_Color C_MAIN     = {255,215,  50, 255};   // gold  — main curve
static constexpr SDL_Color C_OVERLAY  = { 50,200, 255, 255};   // cyan  — second curve
static constexpr SDL_Color C_INTRINS  = { 95,105, 125, 255};   // dim slate — intrinsic
static constexpr SDL_Color C_STRIKE   = {170, 90, 255, 200};   // purple — K line
static constexpr SDL_Color C_BARRIER  = {255,  70,  70, 230};  // red   — H line
static constexpr SDL_Color C_WHITE    = {230, 235, 245, 255};  // near-white text
static constexpr SDL_Color C_DIM      = {105, 118, 140, 255};  // dim text
static constexpr SDL_Color C_GREEN    = { 80, 210, 120, 255};  // green

// =============================================================================
//  Font paths: tried in order; first successful open wins
// =============================================================================
static const char* FONT_PATHS[] = {
    "C:/Windows/Fonts/consola.ttf",          // Consolas (monospace, best for pricer)
    "C:/Windows/Fonts/cour.ttf",             // Courier New
    "C:/Windows/Fonts/arial.ttf",            // Arial
    "C:/Windows/Fonts/calibri.ttf",          // Calibri
    "C:/Windows/Fonts/segoeui.ttf",          // Segoe UI
    "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",   // Linux / MSYS2
    "/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf",
    nullptr
};

TTF_Font* Visualizer::try_load_font(int size) {
    for (int i = 0; FONT_PATHS[i]; ++i) {
        TTF_Font* f = TTF_OpenFont(FONT_PATHS[i], size);
        if (f) return f;
    }
    return nullptr;   // curves still visible without fonts — not fatal
}

// =============================================================================
//  Constructor / Destructor
// =============================================================================
// 32x32 icon: gold call-price curve on dark navy background.
// Generated in pure pixels — no image file dependency.
static SDL_Surface* make_icon() {
    SDL_Surface* s = SDL_CreateRGBSurface(0, 32, 32, 32,
        0x00FF0000, 0x0000FF00, 0x000000FF, 0xFF000000);
    if (!s) return nullptr;

    SDL_FillRect(s, nullptr, SDL_MapRGB(s->format, 13, 17, 23));  // navy bg

    Uint32  gold  = SDL_MapRGB(s->format, 255, 215, 50);
    Uint32* px    = static_cast<Uint32*>(s->pixels);
    int     pitch = s->pitch / 4;

    // Smooth call-price curve: y = 29 - 26 * sigmoid((x-14)/3.5)
    auto cy = [](int x) {
        double sig = 1.0 / (1.0 + std::exp(-(x - 14.0) / 3.5));
        return static_cast<int>(29.0 - 26.0 * sig);
    };

    for (int x = 1; x < 31; ++x) {
        int y0 = std::max(2, std::min(29, cy(x - 1)));
        int y1 = std::max(2, std::min(29, cy(x)));
        for (int y = std::min(y0, y1); y <= std::max(y0, y1); ++y) {
            px[y * pitch + x]     = gold;
            px[y * pitch + x + 1] = gold;   // 2px thick
        }
    }
    return s;
}

Visualizer::Visualizer()  = default;

Visualizer::~Visualizer() {
    if (fontsm_) TTF_CloseFont(fontsm_);
    if (font_)   TTF_CloseFont(font_);
    if (ren_)    SDL_DestroyRenderer(ren_);
    if (win_)    SDL_DestroyWindow(win_);
    if (TTF_WasInit()) TTF_Quit();
    SDL_Quit();
}

// =============================================================================
//  init  — create window + renderer, load fonts
// =============================================================================
bool Visualizer::init() {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        SDL_Log("SDL_Init failed: %s", SDL_GetError());
        return false;
    }

    // TTF not strictly required — continue even if it fails (no text, but curves)
    if (TTF_Init() < 0)
        SDL_Log("SDL_ttf init failed (no text): %s", TTF_GetError());

    win_ = SDL_CreateWindow(
        "BS-Pricer  —  Crank-Nicolson PDE Visualizer",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WIN_W, WIN_H,
        SDL_WINDOW_SHOWN
    );
    if (!win_) { SDL_Log("CreateWindow: %s", SDL_GetError()); return false; }

    ren_ = SDL_CreateRenderer(win_, -1,
              SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!ren_) { SDL_Log("CreateRenderer: %s", SDL_GetError()); return false; }

    SDL_SetRenderDrawBlendMode(ren_, SDL_BLENDMODE_BLEND);

    font_   = try_load_font(16);
    fontsm_ = try_load_font(13);

    SDL_Surface* icon = make_icon();
    if (icon) { SDL_SetWindowIcon(win_, icon); SDL_FreeSurface(icon); }

    return true;
}

// =============================================================================
//  run  — main event loop
// =============================================================================
void Visualizer::run() {
    running_ = true;
    while (running_) {
        handle_events();
        if (dirty_) {
            scene_ = build_scene();   // expensive: runs CN solver
            dirty_ = false;
        }
        render(scene_);   // cheap: always redraw for smooth mouse hover
        SDL_Delay(16);
    }
}

// =============================================================================
//  handle_events  — keyboard input
// =============================================================================
void Visualizer::handle_events() {
    SDL_Event e;
    while (SDL_PollEvent(&e)) {
        if (e.type == SDL_QUIT) { running_ = false; return; }
        if (e.type == SDL_MOUSEMOTION) {
            mouse_x_ = e.motion.x;
            mouse_y_ = e.motion.y;
            continue;
        }
        if (e.type == SDL_MOUSEBUTTONDOWN || e.type == SDL_WINDOWEVENT)
            continue;
        if (e.type != SDL_KEYDOWN) continue;

        bool changed = true;
        switch (e.key.keysym.sym) {
        // --- mode select ---
        case SDLK_1: mode_ = 0; break;
        case SDLK_2: mode_ = 1; break;
        case SDLK_3: mode_ = 2; break;
        case SDLK_4: mode_ = 3; break;
        case SDLK_5: mode_ = 4; break;
        case SDLK_6: mode_ = 5; break;
        // --- sigma ---
        case SDLK_UP:   sigma_ = std::min(sigma_ + 0.01, 0.80); break;
        case SDLK_DOWN: sigma_ = std::max(sigma_ - 0.01, 0.01); break;
        // --- r ---
        case SDLK_RIGHT: r_ = std::min(r_ + 0.005, 0.20);  break;
        case SDLK_LEFT:  r_ = std::max(r_ - 0.005, 0.001); break;
        // --- T ---
        case SDLK_w: T_ = std::min(T_ + 0.1, 5.0); break;
        case SDLK_s: T_ = std::max(T_ - 0.1, 0.1); break;
        // --- strike ---
        case SDLK_d: K_ = std::min(K_ + 5.0, 200.0); break;
        case SDLK_a: K_ = std::max(K_ - 5.0,  10.0); break;
        // --- curve visibility ---
        case SDLK_m: show_main_    = !show_main_;    break;
        case SDLK_o: show_overlay_ = !show_overlay_; break;
        // --- solver toggle ---
        case SDLK_r: use_rcn_ = !use_rcn_; break;
        // --- smoothing toggle (Digital modes only) ---
        case SDLK_g: smooth_ = !smooth_; break;
        // --- quit ---
        case SDLK_ESCAPE: case SDLK_q: running_ = false; changed = false; break;
        default: changed = false; break;
        }
        if (changed) dirty_ = true;
    }
}

// =============================================================================
//  build_scene  — run CN solver(s), assemble CurveData for each mode
// =============================================================================
Visualizer::Scene Visualizer::build_scene() {
    Scene        sc;
    BSParams      p{r_, sigma_, T_};
    CrankNicolson cn_solver;
    ReducedCN     rcn_solver;
    Solver&       solver     = use_rcn_ ? static_cast<Solver&>(rcn_solver)
                                        : static_cast<Solver&>(cn_solver);
    const std::string solver_label = use_rcn_ ? "Reduced CN" : "Crank-Nicolson";

    // Copy solved grid into CurveData (all N+1 nodes)
    auto to_curve = [](const Grid& g) {
        CurveData cd;
        cd.S.reserve(g.N() + 1);
        cd.V.reserve(g.N() + 1);
        for (int i = 0; i <= g.N(); ++i) {
            cd.S.push_back(g.S(i));
            cd.V.push_back(g[i]);
        }
        return cd;
    };

    // Build intrinsic curve (plain max(S-K,0) or max(K-S,0))
    auto make_call_intrinsic = [&](const CurveData& ref) {
        CurveData cd;
        for (double Sv : ref.S) {
            cd.S.push_back(Sv);
            cd.V.push_back(std::max(Sv - K_, 0.0));
        }
        return cd;
    };
    auto make_put_intrinsic = [&](const CurveData& ref) {
        CurveData cd;
        for (double Sv : ref.S) {
            cd.S.push_back(Sv);
            cd.V.push_back(std::max(K_ - Sv, 0.0));
        }
        return cd;
    };

    // === Mode dispatch ==========================================================
    switch (mode_) {

    // ---- [1] European Call  :  FD price + BS analytical overlay ---------------
    case 0: {
        EuropeanCall opt(K_);
        Grid g(SMAX, T_, NV, MV);
        solver.solve(g, opt, p);
        sc.main       = to_curve(g);
        sc.title      = "European Call";
        sc.main_label = "European Call";

        sc.has_overlay   = true;
        sc.overlay_label = "BS Formula";
        for (double Sv : sc.main.S) {
            sc.overlay.S.push_back(Sv);
            sc.overlay.V.push_back(
                Sv > 0.0 ? BSAnalytical::call_price(Sv, K_, r_, sigma_, T_) : 0.0);
        }
        sc.has_intrinsic = true;
        sc.intrinsic = make_call_intrinsic(sc.main);
        break;
    }

    // ---- [2] European Put  :  FD price + BS analytical overlay ----------------
    case 1: {
        EuropeanPut opt(K_);
        Grid g(SMAX, T_, NV, MV);
        solver.solve(g, opt, p);
        sc.main       = to_curve(g);
        sc.title      = "European Put";
        sc.main_label = "European Put";

        sc.has_overlay   = true;
        sc.overlay_label = "BS Formula";
        for (double Sv : sc.main.S) {
            sc.overlay.S.push_back(Sv);
            sc.overlay.V.push_back(BSAnalytical::put_price(Sv, K_, r_, sigma_, T_));
        }
        sc.has_intrinsic = true;
        sc.intrinsic = make_put_intrinsic(sc.main);
        break;
    }

    // ---- [3] Digital Call + Digital Put  :  parity DC + DP = e^{-rT} ---------
    case 2: {
        // smooth_=true  → sigmoid width = dS  (mitigates Gibbs)
        // smooth_=false → width = 1e-6 ≈ step function  (full Gibbs — demo mode)
        const double dS = smooth_ ? (SMAX / NV) : 1e-6;
        DigitalCall dca(K_, dS);
        DigitalPut  dpu(K_, dS);
        Grid gc(SMAX, T_, NV, MV);
        Grid gp(SMAX, T_, NV, MV);
        solver.solve(gc, dca, p);
        solver.solve(gp, dpu, p);
        sc.main          = to_curve(gc);
        sc.title         = smooth_
                           ? "Digital Options  [smoothed - G: show Gibbs]"
                           : "Digital Options  [RAW PAYOFF - Gibbs phenomenon!]";
        sc.main_label    = "Digital Call";
        sc.has_overlay   = true;
        sc.overlay_label = "Digital Put";
        sc.overlay       = to_curve(gp);
        break;
    }

    // ---- [4] American Put  vs  European Put  ----------------------------------
    case 3: {
        AmericanPut am(K_);
        EuropeanPut eu(K_);
        Grid gam(SMAX, T_, NV, MV);
        Grid geu(SMAX, T_, NV, MV);
        solver.solve(gam, am, p);
        solver.solve(geu, eu, p);
        sc.main          = to_curve(gam);
        sc.title         = "American Put vs European Put";
        sc.main_label    = "American Put";
        sc.has_overlay   = true;
        sc.overlay_label = "European Put";
        sc.overlay       = to_curve(geu);
        sc.has_intrinsic = true;
        sc.intrinsic     = make_put_intrinsic(sc.main);
        break;
    }

    // ---- [5] Up-and-Out Call  vs  Vanilla Call  (H = 1.5 * K) ----------------
    case 4: {
        const double H_UP = 1.5 * K_;
        UpAndOutCall ko(K_, H_UP);
        EuropeanCall van(K_);
        Grid gko(H_UP, T_, NV, MV);
        Grid gva(SMAX, T_, NV, MV);
        solver.solve(gko, ko, p);
        solver.solve(gva, van, p);
        sc.main          = to_curve(gko);
        sc.title         = "Up-and-Out Call  (H = 1.5 K)";
        sc.main_label    = "UOC";
        sc.has_overlay   = true;
        sc.overlay_label = "Vanilla Call";
        sc.overlay       = to_curve(gva);
        sc.barrier       = H_UP;
        sc.has_intrinsic = true;
        sc.intrinsic     = make_call_intrinsic(sc.main);
        break;
    }

    // ---- [6] Down-and-Out Put  vs  Vanilla Put  (H = 0.8 * K) ----------------
    case 5: {
        const double H_DN = 0.8 * K_;
        DownAndOutPut ko(K_, H_DN);
        EuropeanPut   van(K_);
        Grid gko(SMAX, T_, NV, MV, H_DN);
        Grid gva(SMAX, T_, NV, MV);
        solver.solve(gko, ko, p);
        solver.solve(gva, van, p);
        sc.main          = to_curve(gko);
        sc.title         = "Down-and-Out Put  (H = 0.8 K)";
        sc.main_label    = "DOP";
        sc.has_overlay   = true;
        sc.overlay_label = "Vanilla Put";
        sc.overlay       = to_curve(gva);
        sc.barrier       = H_DN;
        sc.has_intrinsic = true;
        sc.intrinsic     = make_put_intrinsic(sc.main);
        break;
    }
    }

    // === Compute V range ========================================================
    sc.S_lo = 0.0;
    sc.S_hi = SMAX;   // fixed axis so K line visibly shifts when K changes

    double vmax = 5.0;    // minimum scale (avoid degenerate all-zero chart)
    auto update_vmax = [&](const CurveData& cd) {
        for (size_t i = 0; i < cd.S.size(); ++i)
            if (cd.S[i] >= sc.S_lo && cd.S[i] <= sc.S_hi)
                vmax = std::max(vmax, cd.V[i]);
    };
    update_vmax(sc.main);
    if (sc.has_overlay)   update_vmax(sc.overlay);
    if (sc.has_intrinsic) update_vmax(sc.intrinsic);

    sc.V_lo = 0.0;
    sc.V_hi = vmax * 1.15;

    return sc;
}

// =============================================================================
//  Coordinate transforms  (world → screen pixels)
// =============================================================================
int Visualizer::sx(double S, const Scene& sc) const {
    double t = (S - sc.S_lo) / (sc.S_hi - sc.S_lo);
    return PX0 + static_cast<int>(t * PW + 0.5);
}

int Visualizer::vx(double V, const Scene& sc) const {
    double t = (V - sc.V_lo) / (sc.V_hi - sc.V_lo);
    int y    = PY0 + PH - static_cast<int>(t * PH + 0.5);
    return std::max(PY0 - 5, std::min(PY0 + PH + 5, y));   // small overdraw allowed
}

// =============================================================================
//  Drawing primitives
// =============================================================================
void Visualizer::fill_rect(int x, int y, int w, int h, SDL_Color c) {
    SDL_SetRenderDrawColor(ren_, c.r, c.g, c.b, c.a);
    SDL_Rect r{x, y, w, h};
    SDL_RenderFillRect(ren_, &r);
}

void Visualizer::draw_line(int x1, int y1, int x2, int y2, SDL_Color c) {
    SDL_SetRenderDrawColor(ren_, c.r, c.g, c.b, c.a);
    SDL_RenderDrawLine(ren_, x1, y1, x2, y2);
}

void Visualizer::draw_text(const std::string& s, int x, int y,
                            SDL_Color col, bool large) {
    TTF_Font* f = large ? font_ : fontsm_;
    if (!f || s.empty()) return;
    SDL_Surface* surf = TTF_RenderText_Blended(f, s.c_str(), col);
    if (!surf) return;
    SDL_Texture* tex = SDL_CreateTextureFromSurface(ren_, surf);
    if (tex) {
        SDL_Rect dst{x, y, surf->w, surf->h};
        SDL_RenderCopy(ren_, tex, nullptr, &dst);
        SDL_DestroyTexture(tex);
    }
    SDL_FreeSurface(surf);
}

// Draw a curve; thickness = 1 or 2 (repeat-draw offset for thick lines).
// dashed = skip every other segment (for intrinsic / reference lines).
void Visualizer::draw_curve(const CurveData& cd, const Scene& sc,
                              SDL_Color col, int thickness, bool dashed) {
    SDL_SetRenderDrawColor(ren_, col.r, col.g, col.b, col.a);
    int prev_x = -1, prev_y = -1;
    int seg = 0;

    for (size_t i = 0; i < cd.S.size(); ++i) {
        double Sv = cd.S[i];
        if (Sv < sc.S_lo || Sv > sc.S_hi) { prev_x = -1; continue; }
        int cx = sx(Sv, sc);
        int cy = vx(cd.V[i], sc);
        if (prev_x >= 0) {
            if (!dashed || (seg % 6) < 4) {   // 4-on, 2-off for dashed
                SDL_RenderDrawLine(ren_, prev_x, prev_y, cx, cy);
                if (thickness >= 2) {
                    SDL_RenderDrawLine(ren_, prev_x, prev_y + 1, cx, cy + 1);
                }
            }
            ++seg;
        }
        prev_x = cx;
        prev_y = cy;
    }
}

// Dashed vertical line at world-space S.
void Visualizer::draw_vline_dashed(double S, const Scene& sc, SDL_Color col) {
    if (S < sc.S_lo || S > sc.S_hi) return;
    int x = sx(S, sc);
    SDL_SetRenderDrawColor(ren_, col.r, col.g, col.b, col.a);
    for (int y = PY0; y < PY0 + PH; y += 9) {
        int y2 = std::min(y + 5, PY0 + PH);
        SDL_RenderDrawLine(ren_, x,     y,  x,     y2);
        SDL_RenderDrawLine(ren_, x + 1, y,  x + 1, y2);
    }
}

// =============================================================================
//  draw_axes  —  background, grid lines, axis tick labels
// =============================================================================
void Visualizer::draw_axes(const Scene& sc) {
    // Plot background
    fill_rect(PX0, PY0, PW, PH, C_PLOTBG);

    // Horizontal grid lines (V axis)
    constexpr int N_HGRID = 5;
    for (int k = 0; k <= N_HGRID; ++k) {
        double V = sc.V_lo + k * (sc.V_hi - sc.V_lo) / N_HGRID;
        int    y = vx(V, sc);
        SDL_SetRenderDrawColor(ren_, C_GRIDLN.r, C_GRIDLN.g, C_GRIDLN.b, 255);
        SDL_RenderDrawLine(ren_, PX0, y, PX0 + PW, y);
        // Y label
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << V;
        draw_text(oss.str(), PX0 - 52, y - 7, C_DIM);
    }

    // Vertical grid lines (S axis)
    constexpr int N_VGRID = 8;
    for (int k = 0; k <= N_VGRID; ++k) {
        double S = sc.S_lo + k * (sc.S_hi - sc.S_lo) / N_VGRID;
        int    x = sx(S, sc);
        SDL_SetRenderDrawColor(ren_, C_GRIDLN.r, C_GRIDLN.g, C_GRIDLN.b, 255);
        SDL_RenderDrawLine(ren_, x, PY0, x, PY0 + PH);
        // X label
        std::ostringstream oss;
        oss << static_cast<int>(std::round(S));
        draw_text(oss.str(), x - 12, PY0 + PH + 6, C_DIM);
    }

    // Axis border (L + bottom)
    SDL_SetRenderDrawColor(ren_, C_DIM.r, C_DIM.g, C_DIM.b, 255);
    SDL_RenderDrawLine(ren_, PX0, PY0, PX0, PY0 + PH);
    SDL_RenderDrawLine(ren_, PX0, PY0 + PH, PX0 + PW, PY0 + PH);

    // Axis labels
    draw_text("S", PX0 + PW / 2 - 4, PY0 + PH + 22, C_DIM);
    draw_text("V", PX0 - 65, PY0 + PH / 2, C_DIM);
}

// =============================================================================
//  draw_side_panel  —  right-hand parameters + ATM price panel
// =============================================================================
void Visualizer::draw_side_panel(const Scene& sc) {
    const int xp = PX0 + PW + 22;
    const int wp = WIN_W - xp - 10;

    fill_rect(xp - 8, PY0, wp + 8, WIN_H - PY0 - 10, C_PANEL);

    int yp = PY0 + 14;
    constexpr int LH = 23;   // line height

    draw_text("PARAMETERS", xp, yp, C_WHITE, true); yp += LH + 4;

    // separator
    SDL_SetRenderDrawColor(ren_, 45, 60, 90, 255);
    SDL_RenderDrawLine(ren_, xp, yp, xp + wp - 16, yp);
    yp += 10;

    // Helper to print "label   value  [key]"
    auto param_row = [&](const std::string& lbl, double val, int dec,
                          const std::string& key, SDL_Color vcol) {
        draw_text(lbl, xp, yp, C_DIM);
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(dec) << val;
        draw_text(oss.str(), xp + 65, yp, vcol);
        draw_text(key, xp + 130, yp, C_DIM);
        yp += LH;
    };

    param_row("sigma:", sigma_, 3, "[up/dn]", C_MAIN);
    param_row("r:",     r_,     3, "[lt/rt]", C_MAIN);
    param_row("T:",     T_,     2, "[W/S]",   C_MAIN);
    param_row("K:",     K_,     0, "[A/D]",   C_MAIN);

    // Solver name
    draw_text("solver:", xp, yp, C_DIM);
    draw_text(use_rcn_ ? "Reduced CN" : "Crank-Nicols.", xp + 65, yp, C_GREEN);
    draw_text("[R]", xp + 185, yp, C_DIM);
    yp += LH;

    yp += 6;
    SDL_SetRenderDrawColor(ren_, 45, 60, 90, 255);
    SDL_RenderDrawLine(ren_, xp, yp, xp + wp - 16, yp);
    yp += 10;

    draw_text("ATM  (S = K)", xp, yp, C_WHITE, true); yp += LH + 4;

    // Interpolate a CurveData at S = K_
    auto interp_at_K = [&](const CurveData& cd) -> double {
        for (size_t i = 0; i + 1 < cd.S.size(); ++i) {
            if (cd.S[i] <= K_ && cd.S[i + 1] >= K_) {
                double a = (K_ - cd.S[i]) / (cd.S[i + 1] - cd.S[i]);
                return (1.0 - a) * cd.V[i] + a * cd.V[i + 1];
            }
        }
        return std::numeric_limits<double>::quiet_NaN();
    };

    auto print_atm = [&](const std::string& lbl, const CurveData& cd, SDL_Color col) {
        double v = interp_at_K(cd);
        std::ostringstream oss;
        if (std::isnan(v)) oss << "-";
        else               oss << std::fixed << std::setprecision(4) << v;
        draw_text(lbl + ":", xp, yp, C_DIM);
        draw_text(oss.str(), xp + 118, yp, col);
        yp += LH + 2;
    };

    // Modes 0/1: solver name is the meaningful label (CN vs ReducedCN comparison)
    // Other modes: the option type is more descriptive
    print_atm(sc.main_label, sc.main, C_MAIN);
    if (sc.has_overlay)
        print_atm(sc.overlay_label, sc.overlay, C_OVERLAY);

    // For European Call/Put show the BS ATM price explicitly
    auto stat_row = [&](const std::string& lbl, const std::string& val) {
        draw_text(lbl, xp, yp, C_DIM);
        draw_text(val, xp + 118, yp, C_DIM);
        yp += LH + 2;
    };

    if (mode_ == 0) {
        double bs_atm = BSAnalytical::call_price(K_, K_, r_, sigma_, T_);
        double fd_atm = interp_at_K(sc.main);
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(2) << std::abs(fd_atm - bs_atm);
        stat_row("Err FD-BS:", oss.str());
    } else if (mode_ == 1) {
        double bs_atm = BSAnalytical::put_price(K_, K_, r_, sigma_, T_);
        double fd_atm = interp_at_K(sc.main);
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(2) << std::abs(fd_atm - bs_atm);
        stat_row("Err FD-BS:", oss.str());
    } else if (mode_ == 2) {
        double dc = interp_at_K(sc.main);
        double dp = interp_at_K(sc.overlay);
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(5) << (dc + dp);
        stat_row("DC+DP:", oss.str());
        std::ostringstream oss2;
        oss2 << std::fixed << std::setprecision(5) << std::exp(-r_ * T_);
        stat_row("e^-rT:", oss2.str());
    } else if (mode_ == 3) {
        double am = interp_at_K(sc.main);
        double eu = interp_at_K(sc.overlay);
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4) << (am - eu);
        stat_row("EE Premium:", oss.str());
    }

    // Mode list
    yp += 10;
    SDL_SetRenderDrawColor(ren_, 45, 60, 90, 255);
    SDL_RenderDrawLine(ren_, xp, yp, xp + wp - 16, yp);
    yp += 10;
    draw_text("MODES", xp, yp, C_WHITE, true); yp += LH + 2;

    static const char* MODE_NAMES[] = {
        "1  EU Call",
        "2  EU Put",
        "3  Digital",
        "4  AM Put",
        "5  UO Call",
        "6  DO Put",
    };
    for (int m = 0; m < 6; ++m) {
        SDL_Color col = (m == mode_) ? C_MAIN : C_DIM;
        draw_text(MODE_NAMES[m], xp, yp, col);
        yp += LH - 2;
    }
}

// =============================================================================
//  draw_legend  —  below the plot area, above the help bar
// =============================================================================
void Visualizer::draw_legend(const Scene& sc) {
    // Position: just below the X axis tick labels
    const int lx = PX0;
    const int ly = PY0 + PH + 34;   // below tick labels (~22px) + small gap

    int x = lx;
    auto entry = [&](SDL_Color col, const std::string& label, bool dashed) {
        const int MID = ly + 7;
        if (dashed) {
            SDL_SetRenderDrawColor(ren_, col.r, col.g, col.b, col.a);
            for (int dx = 0; dx < 16; dx += 5)
                SDL_RenderDrawLine(ren_, x + dx, MID, x + dx + 3, MID);
        } else {
            SDL_SetRenderDrawColor(ren_, col.r, col.g, col.b, col.a);
            for (int dy = MID - 2; dy <= MID + 2; ++dy)
                SDL_RenderDrawLine(ren_, x, dy, x + 16, dy);
        }
        draw_text(label, x + 20, ly, col);
        x += 20 + static_cast<int>(label.size()) * 7 + 18;
    };

    entry(show_main_    ? C_MAIN    : C_DIM, sc.main_label, false);
    if (sc.has_overlay)
        entry(show_overlay_ ? C_OVERLAY : C_DIM, sc.overlay_label, false);
    if (sc.has_intrinsic)
        entry(C_INTRINS, "Intrinsic", true);
    if (sc.barrier > 0)
        entry(C_BARRIER, "Barrier H", true);
}

// =============================================================================
//  draw_help  —  bottom keyboard shortcut bar
// =============================================================================
// =============================================================================
//  draw_hover  —  crosshair + tooltip when mouse is inside the plot area
// =============================================================================
void Visualizer::draw_hover(const Scene& sc) {
    // Only active inside the plot rectangle, and only when at least one curve is visible
    if (mouse_x_ < PX0 || mouse_x_ >= PX0 + PW ||
        mouse_y_ < PY0 || mouse_y_ >= PY0 + PH)
        return;
    if (!show_main_ && !show_overlay_)
        return;

    // World coordinates under the cursor
    double S_h = sc.S_lo + (mouse_x_ - PX0) * (sc.S_hi - sc.S_lo) / PW;

    // Interpolate the main curve at S_h
    auto interp = [&](const CurveData& cd) -> double {
        for (size_t i = 0; i + 1 < cd.S.size(); ++i) {
            if (cd.S[i] <= S_h && cd.S[i + 1] >= S_h) {
                double a = (S_h - cd.S[i]) / (cd.S[i + 1] - cd.S[i]);
                return (1.0 - a) * cd.V[i] + a * cd.V[i + 1];
            }
        }
        return std::numeric_limits<double>::quiet_NaN();
    };

    double V_main = interp(sc.main);
    double V_over = sc.has_overlay ? interp(sc.overlay)
                                   : std::numeric_limits<double>::quiet_NaN();

    // --- Crosshair lines (clipped to plot area) ---
    SDL_Rect clip{PX0, PY0, PW, PH};
    SDL_RenderSetClipRect(ren_, &clip);

    SDL_SetRenderDrawColor(ren_, 160, 165, 185, 90);
    SDL_RenderDrawLine(ren_, mouse_x_, PY0, mouse_x_, PY0 + PH);   // vertical
    SDL_RenderDrawLine(ren_, PX0, mouse_y_, PX0 + PW, mouse_y_);   // horizontal

    // Dot on the main curve
    if (show_main_ && !std::isnan(V_main)) {
        int cy = vx(V_main, sc);
        SDL_SetRenderDrawColor(ren_, C_MAIN.r, C_MAIN.g, C_MAIN.b, 255);
        for (int dx = -3; dx <= 3; ++dx)
            for (int dy = -3; dy <= 3; ++dy)
                if (dx*dx + dy*dy <= 9)
                    SDL_RenderDrawPoint(ren_, mouse_x_ + dx, cy + dy);
    }
    // Dot on the overlay curve
    if (!std::isnan(V_over) && show_overlay_) {
        int cy = vx(V_over, sc);
        SDL_SetRenderDrawColor(ren_, C_OVERLAY.r, C_OVERLAY.g, C_OVERLAY.b, 255);
        for (int dx = -3; dx <= 3; ++dx)
            for (int dy = -3; dy <= 3; ++dy)
                if (dx*dx + dy*dy <= 9)
                    SDL_RenderDrawPoint(ren_, mouse_x_ + dx, cy + dy);
    }

    SDL_RenderSetClipRect(ren_, nullptr);

    // --- Tooltip box ---
    std::ostringstream line1, line2, line3;
    line1 << "S = " << std::fixed << std::setprecision(2) << S_h;
    if (show_main_ && !std::isnan(V_main))
        line2 << sc.main_label << " = " << std::fixed << std::setprecision(4) << V_main;
    if (!std::isnan(V_over) && show_overlay_)
        line3 << sc.overlay_label << " = " << std::fixed << std::setprecision(4) << V_over;

    const int TW = 165, TH = line3.str().empty() ? 38 : 56;
    int tx = mouse_x_ + 12;
    int ty = mouse_y_ - TH - 6;
    if (tx + TW > PX0 + PW) tx = mouse_x_ - TW - 8;
    if (ty < PY0 + 4)       ty = mouse_y_ + 10;

    fill_rect(tx - 4, ty - 3, TW, TH, {13, 17, 30, 220});
    SDL_SetRenderDrawColor(ren_, 60, 80, 120, 200);
    SDL_Rect border{tx - 4, ty - 3, TW, TH};
    SDL_RenderDrawRect(ren_, &border);

    draw_text(line1.str(), tx, ty,      C_WHITE);
    if (!line2.str().empty()) draw_text(line2.str(), tx, ty + 16, C_MAIN);
    if (!line3.str().empty()) draw_text(line3.str(), tx, ty + 32, C_OVERLAY);
}

void Visualizer::draw_help() {
    const int y = WIN_H - 20;
    SDL_SetRenderDrawColor(ren_, 25, 30, 50, 255);
    SDL_Rect bar{0, y - 4, WIN_W, 24};
    SDL_RenderFillRect(ren_, &bar);
    std::string help =
        "1-6: mode   up/dn: sigma   lt/rt: r   W/S: T   A/D: K   R: solver   M/O: curves";
    if (mode_ == 2)
        help += std::string("   G: ") + (smooth_ ? "show Gibbs" : "smooth");
    help += "   ESC: quit";
    draw_text(help, PX0, y, C_DIM);
}

// =============================================================================
//  render  —  compose full frame from the pre-built Scene
// =============================================================================
void Visualizer::render(const Scene& sc) {
    // Background
    SDL_SetRenderDrawColor(ren_, C_BG.r, C_BG.g, C_BG.b, 255);
    SDL_RenderClear(ren_);

    // Axes + grid lines
    draw_axes(sc);

    // Strike line  (K, dim yellow dashed)
    // Vertical reference lines (drawn before clip so they show fully)
    draw_vline_dashed(K_, sc, C_STRIKE);
    if (sc.barrier > 0.0)
        draw_vline_dashed(sc.barrier, sc, C_BARRIER);

    // ---- Clip curves to plot area (SDL native clip — does NOT cover axis labels) ----
    SDL_Rect clip_rect{PX0, PY0, PW, PH};
    SDL_RenderSetClipRect(ren_, &clip_rect);

    if (sc.has_intrinsic)
        draw_curve(sc.intrinsic, sc, C_INTRINS, 1, true);
    if (sc.has_overlay && show_overlay_)
        draw_curve(sc.overlay,   sc, C_OVERLAY, 2, false);
    if (show_main_)
        draw_curve(sc.main,      sc, C_MAIN,    2, false);

    SDL_RenderSetClipRect(ren_, nullptr);   // remove clip
    // ---- End clip ---------------------------------------------------------------

    // Axis border (drawn after curves so it sits on top)
    SDL_SetRenderDrawColor(ren_, C_DIM.r, C_DIM.g, C_DIM.b, 255);
    SDL_RenderDrawLine(ren_, PX0,      PY0,      PX0,      PY0 + PH);
    SDL_RenderDrawLine(ren_, PX0,      PY0 + PH, PX0 + PW, PY0 + PH);

    // K / H labels (inside plot, above axis line)
    {
        int kx = sx(K_, sc);
        std::ostringstream oss;
        oss << "K=" << static_cast<int>(K_);
        draw_text(oss.str(), kx + 4, PY0 + 4, C_STRIKE);
    }
    if (sc.barrier > 0.0) {
        int bx = sx(sc.barrier, sc);
        std::ostringstream oss;
        oss << "H=" << static_cast<int>(std::round(sc.barrier));
        draw_text(oss.str(), bx + 4, PY0 + 18, C_BARRIER);
    }

    // Legend — below the plot, above the help bar
    draw_legend(sc);

    // Title
    {
        std::string title = "[" + std::to_string(mode_ + 1) + "]  " + sc.title;
        SDL_Color title_col = (mode_ == 2 && !smooth_) ? C_BARRIER : C_WHITE;
        draw_text(title, PX0, 16, title_col, true);
    }

    // Side panel
    draw_side_panel(sc);

    // Hover crosshair + tooltip
    draw_hover(sc);

    // Help bar
    draw_help();

    SDL_RenderPresent(ren_);
}
