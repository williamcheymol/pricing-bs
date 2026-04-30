CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra

# ---- pricer (console, no SDL) -----------------------------------------------
SRCS := Grid.cpp CrankNicolson.cpp ReducedCN.cpp Greeks.cpp BSFormula.cpp Thomas.cpp main.cpp
OBJS := $(SRCS:.cpp=.o)

.PHONY: all viz clean

all: pricer

pricer: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---- pricer_viz (SDL2 interactive window) ------------------------------------
#
#  Requires:  SDL2  +  SDL2_ttf
#  MSYS2/MinGW:  pacman -S mingw-w64-x86_64-SDL2 mingw-w64-x86_64-SDL2_ttf
#  Ubuntu/Debian: apt  install libsdl2-dev libsdl2-ttf-dev
#  macOS (brew):  brew install sdl2 sdl2_ttf
#
#  sdl2-config is the canonical way to get compiler / linker flags; fall back
#  to empty strings if it is not on PATH (compilation will still fail loudly
#  with a helpful missing-header error rather than a cryptic linker error).

# sdl2-config is in the MSYS2 shell PATH but not always in the Windows PATH.
# We try it first; if it returns nothing we fall back to the known MSYS2 paths.
_SDL_CF := $(shell C:/msys64/mingw64/bin/sdl2-config --cflags 2>/dev/null)
_SDL_LF := $(shell C:/msys64/mingw64/bin/sdl2-config --libs   2>/dev/null)

SDL_CFLAGS := $(if $(_SDL_CF),$(_SDL_CF),-IC:/msys64/mingw64/include)
SDL_LIBS   := $(if $(_SDL_LF),$(_SDL_LF),-LC:/msys64/mingw64/lib -lmingw32 -lSDL2main -lSDL2) -lSDL2_ttf

# SDL2 was installed into the MSYS2/MinGW-w64 toolchain.
# Use the MSYS2 64-bit g++ explicitly so headers and libs match.
VIZ_CXX  := C:/msys64/mingw64/bin/g++
VIZ_FLAGS := -std=c++17 -O2 -Wall -Wextra

VIZ_SHARED := Grid.cpp CrankNicolson.cpp ReducedCN.cpp Greeks.cpp \
              BSFormula.cpp Thomas.cpp
VIZ_SRCS   := $(VIZ_SHARED) Visualizer.cpp main_viz.cpp

viz: $(VIZ_SRCS)
	$(VIZ_CXX) $(VIZ_FLAGS) $(SDL_CFLAGS) $^ $(SDL_LIBS) -o pricer_viz

# ---- clean -------------------------------------------------------------------
clean:
	rm -f $(OBJS) pricer pricer.exe pricer_viz pricer_viz.exe
