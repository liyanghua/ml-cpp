CC = gcc
CXX = g++
AR = ar
RANLIB = ranlib
CPPFLAGS =
CXXFLAGS = -Wall -g -O3
LIBS =

all: owlqn

owlqn: leastSquares.o logreg.o main.o OWLQN.o TerminationCriterion.o
	$(CXX) -o $@ $?

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.a owlqn
