# compiler settings
CC = g++
CFLAGS = -std=c++17 -fopenmp -Wpedantic -Wall -O3

SRC_FILES = $(wildcard *.cpp)
OBJ_FILES = $(patsubst %.cpp, %.o ,$(SRC_FILES))

all : UncertainCentrality
.PHONY : all

UncertainCentrality: $(OBJ_FILES)
	$(CC) $(CFLAGS) -o UncertainCentrality $^

main.o : main.cpp
	$(CC) $(CFLAGS) -c -o main.o main.cpp

UncertainGraph.o : UncertainGraph.cpp
	$(CC) $(CFLAGS)  -c -o UncertainGraph.o UncertainGraph.cpp

UncertainGraphAlgorithms.o : UncertainGraphAlgorithms.cpp
	$(CC) $(CFLAGS)  -c -o UncertainGraphAlgorithms.o UncertainGraphAlgorithms.cpp


clean :
	rm *.o
.PHONY : clean
