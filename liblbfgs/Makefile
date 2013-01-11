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
	$(AR) -rc $@ $?
	$(RANLIB) $@

sample_c: sample_c.o liblbfgs.a
	$(CC) $(LIBS) -o $@ $?

sample_cpp: sample_cpp.o liblbfgs.a
	$(CXX) $(LIBS) -o $@ $?

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.a sample_c sample_cpp
