CC = gcc
CXX = g++
AR = ar
RANLIB = ranlib
CPPFLAGS = -Isrc
CXXFLAGS = -g -O3
LIBS =

all: libalglib.a test_c test_i

libalglib.a: src/alglibinternal.o src/alglibmisc.o src/ap.o src/dataanalysis.o src/diffequations.o src/fasttransforms.o src/integration.o src/interpolation.o src/linalg.o src/optimization.o src/solvers.o src/specialfunctions.o src/statistics.o
	$(AR) -rc $@ $?
	$(RANLIB) $@

test_c: tests/test_c.o libalglib.a
	$(CXX) $(LIBS) -o $@ $?

test_i: tests/test_i.o libalglib.a
	$(CXX) $(LIBS) -o $@ $?

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f src/*.o tests/*.o *.a test_c test_i
