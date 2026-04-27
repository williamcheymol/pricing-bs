CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra

SRCS := Grid.cpp CrankNicolson.cpp ReducedCN.cpp Greeks.cpp BSFormula.cpp Thomas.cpp main.cpp
OBJS := $(SRCS:.cpp=.o)

.PHONY: all clean

all: pricer

pricer: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) pricer pricer.exe
