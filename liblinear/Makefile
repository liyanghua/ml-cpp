CC = gcc
CXX = g++
AR = ar
RANLIB = ranlib
CPPFLAGS =
CFLAGS = -Wall -g -O3
CXXFLAGS = $(CFLAGS)
LIBS =

all: liblinear.a train predict

liblinear.a: blas/daxpy.o blas/ddot.o blas/dnrm2.o blas/dscal.o linear.o tron.o
	$(AR) -rc $@ $?
	$(RANLIB) $@

train: train.o liblinear.a
	$(CXX) $(LIBS) -o $@ $?

predict: predict.o liblinear.a
	$(CXX) $(LIBS) -o $@ $?

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.a train predict
