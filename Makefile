CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -Isrc

# ---- pricer (console, no SDL) -----------------------------------------------
SRCS := src/Grid.cpp src/CrankNicolson.cpp src/ReducedCN.cpp src/Greeks.cpp \
        src/BSFormula.cpp src/Thomas.cpp main.cpp
OBJS := $(notdir $(SRCS:.cpp=.o))

.PHONY: all viz clean

all: pricer

pricer: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---- pricer_viz (SDL2 interactive window) ------------------------------------
#
#  Requires:  SDL2  +  SDL2_ttf
#  MSYS2/MinGW:  pacman -S mingw-w64-x86_64-SDL2 mingw-w64-x86_64-SDL2_ttf
#  Ubuntu/Debian: apt  install libsdl2-dev libsdl2-ttf-dev
#  macOS (brew):  brew install sdl2 sdl2_ttf

_SDL_CF := $(shell C:/msys64/mingw64/bin/sdl2-config --cflags 2>/dev/null)
_SDL_LF := $(shell C:/msys64/mingw64/bin/sdl2-config --libs   2>/dev/null)

SDL_CFLAGS := $(if $(_SDL_CF),$(_SDL_CF),-IC:/msys64/mingw64/include)
SDL_LIBS   := $(if $(_SDL_LF),$(_SDL_LF),-LC:/msys64/mingw64/lib -lmingw32 -lSDL2main -lSDL2) -lSDL2_ttf

VIZ_CXX   := C:/msys64/mingw64/bin/g++
VIZ_FLAGS := -std=c++17 -O2 -Wall -Wextra -Isrc -Iviz

VIZ_SRCS := src/Grid.cpp src/CrankNicolson.cpp src/ReducedCN.cpp src/Greeks.cpp \
            src/BSFormula.cpp src/Thomas.cpp viz/Visualizer.cpp viz/main_viz.cpp

viz: $(VIZ_SRCS)
	$(VIZ_CXX) $(VIZ_FLAGS) $(SDL_CFLAGS) $^ $(SDL_LIBS) -o pricer_viz

# ---- clean -------------------------------------------------------------------
clean:
	rm -f $(OBJS) pricer pricer.exe pricer_viz pricer_viz.exe
