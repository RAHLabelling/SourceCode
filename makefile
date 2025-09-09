CC = g++ -std=c++20 -O3 -Wall -Wextra -o
TCC = g++ -std=c++20 -ggdb -Wall -Wextra -o
INC = road_network.cpp util.cpp

default: main
all: main index query topcut test
main:
	$(CC) cut main.cpp $(INC)
index:
	$(CC) index index.cpp $(INC)
query:
	$(CC) query query.cpp $(INC)
topcut:
	$(CC) topcut topcut.cpp $(INC)
test:
	$(TCC) test test.cpp $(INC)
clean:
	rm cut index query topcut test

.PHONY: main index query topcut test
