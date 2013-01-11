CC = gcc
CXX = g++
AR = ar
RANLIB = ranlib
CPPFLAGS = -DUSE_SSE=1 -DHAVE_EMMINTRIN_H=1
CFLAGS = -Wall -g -O3
CXXFLAGS = $(CFLAGS)
LIBS = -lm

all: liblbfgs.a sample_c sample_cpp

liblbfgs.a: lbfgs.o
	$(AR) -rc $@ lbfgs.o
	$(RANLIB) $@

sample_c: sample_c.o
	$(CC) $(LIBS) -o $@ $< liblbfgs.a

sample_cpp: sample_cpp.o
	$(CXX) $(LIBS) -o $@ $< liblbfgs.a

sample_c.o: sample.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

sample_cpp.o: sample.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.a sample_c sample_cpp
