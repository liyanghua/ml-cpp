include ../Makefile.in

all: example1 example2

example1: 
	$(CC) -c $(CFLAGS) -I../include -I../src cblas_example1.c
	$(LOADER) -o cblas_ex1 cblas_example1.o $(CBLIB) $(BLLIB)

example2: 
	$(CC) -c $(CFLAGS) -I../include -I../src cblas_example2.c
	$(LOADER) -o cblas_ex2 cblas_example2.o $(CBLIB) $(BLLIB)
   
clean:
	rm -f *.o cblas_ex1 cblas_ex2
