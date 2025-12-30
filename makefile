# Compiler
CXX := g++

# Release / Debug flags
CXXFLAGS  := -std=c++20 -O3 -Wall -Wextra -fopenmp -mavx2 -mfma -pthread
DBGFLAGS  := -std=c++20 -ggdb -O0 -Wall -Wextra -fopenmp -mavx2 -mfma -pthread

CXXFLAGS  += -flto
DBGFLAGS  += -flto

# Sources
INC := road_network.cpp util.cpp

default: main
all: main index query topcut test

main: main.cpp $(INC)
	$(CXX) $(CXXFLAGS) -o cut main.cpp $(INC)

index: index.cpp $(INC)
	$(CXX) $(CXXFLAGS) -o index index.cpp $(INC)

query: query.cpp $(INC)
	$(CXX) $(CXXFLAGS) -o query query.cpp $(INC)

topcut: topcut.cpp $(INC)
	$(CXX) $(CXXFLAGS) -o topcut topcut.cpp $(INC)

test: test.cpp $(INC)
	$(CXX) $(DBGFLAGS) -o test test.cpp $(INC)

clean:
	rm -f cut index query topcut test

.PHONY: default all main index query topcut test clean
