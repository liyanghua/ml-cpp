dlvl = ./.
include $(dlvl)/Makefile.in

all: alllib alltst

clean:
	( cd testing && make clean )
	( cd src && make clean )
	( cd examples && make clean )

slib1:  sreal1
dlib1:  dreal1
clib1:  scplx1
zlib1:  dcplx1
slib2:  sreal2
dlib2:  dreal2
clib2:  scplx2
zlib2:  dcplx2
slib3:  sreal3
dlib3:  dreal3
clib3:  scplx3
zlib3:  dcplx3
alllib1: allprecision1
alllib2: allprecision2
alllib3: allprecision3
alllib:  allprecision


sreal1:
	( cd src && make slib1)
dreal1:
	( cd src && make dlib1)
scplx1:
	( cd src && make clib1)
dcplx1:
	( cd src && make zlib1)
allprecision1:
	( cd src && make all1)
sreal2:
	( cd src && make slib2)
dreal2:
	( cd src && make dlib2)
scplx2:
	( cd src && make clib2)
dcplx2:
	( cd src && make zlib2)
allprecision2:
	( cd src && make all2)
sreal3:
	( cd src && make slib3)
dreal3:
	( cd src && make dlib3)
scplx3:
	( cd src && make clib3)
dcplx3:
	( cd src && make zlib3)
allprecision3:
	( cd src && make all3)
allprecision:
	( cd src && make all)

stest1:
	( cd testing && make stest1 )
dtest1:
	( cd testing && make dtest1 )
ctest1:
	( cd testing && make ctest1 )
ztest1:
	( cd testing && make ztest1 )
alltst1:
	( cd testing && make all1 )
stest2:
	( cd testing && make stest2 )
dtest2:
	( cd testing && make dtest2 )
ctest2:
	( cd testing && make ctest2 )
ztest2:
	( cd testing && make ztest2 )
alltst2:
	( cd testing && make all2 )
stest3:
	( cd testing && make stest3 )
dtest3:
	( cd testing && make dtest3 )
ctest3:
	( cd testing && make ctest3 )
ztest3:
	( cd testing && make ztest3 )
alltst3:
	( cd testing && make all3 )
alltst:
	( cd testing && make all )
runtst:
	( cd testing && make run )

example: alllib
	( cd examples && make all )

