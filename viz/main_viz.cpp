#include "Visualizer.hpp"

// =============================================================================
//  main_viz.cpp  —  Entry point for the SDL2 interactive visualizer
//
//  Build:   make viz
//  Run:     ./pricer_viz    (or pricer_viz.exe on Windows)
//
//  Controls:
//    1 – 6    switch display mode
//    ↑ / ↓   sigma  ± 0.01
//    ← / →   r      ± 0.005
//    W / S    T      ± 0.1 years
//    ESC/Q    quit
// =============================================================================

// SDL2 on Windows redefines main as SDL_main via the macro in SDL.h.
// Providing argc/argv keeps the linker happy with SDL2main.
int main(int /*argc*/, char* /*argv*/[]) {
    Visualizer viz;
    if (!viz.init()) return 1;
    viz.run();
    return 0;
}
