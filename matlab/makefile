# This is the makefile to create the mex-function for lu1fac
# from file lu1fac.c, mi27lu.f, mi15blas.f and sort.c.
#
# 30 Apr 2004: (or earlier) Original version by Michael O'Sullivan
# 18 Jan 2008: Switched mex to gfortran

# The matlab directory
MATLAB =
# The matlab include file directory
INCLUDE = /usr/include

# The libraries
LIBS =

# The compile options
#OPTIONS  = -g -O
#F77OPTS  = -g -c -O -64 -i4
OPTIONS  = -O -v
#F77OPTS  = -c -O -64 -i4
F77OPTS  = -c -O -i4

# Clear suffix list, then define the ones we want
.SUFFIXES:
.SUFFIXES: .c .f .f90 .o

.f90.o:; ${FC} ${FFLAGS} -c -o $@ $<
.f.o:;   ${FC} ${FFLAGS} -c -o $@ $<
.c.o:;   $(CC) $(CFLAGS) -c -o $@ $<

#lu1fac.mexsg : lu1fac.c mi27lu.o mi15blas.o sort.c
#	mex $(OPTIONS)  -L/usr/lib/gcc-lib/i486-suse-linux/2.95.3 -lg2c lu1fac.c -I$(INCLUDE) mi27lu.o mi15blas.o sort.c $(LIBS)

#lu1fac.mexglx : lu1fac.c mi27lu.o mi15blas.o sort.c
#	mex $(OPTIONS)  -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2 -lg2c lu1fac.c -I$(INCLUDE) mi27lu.o mi15blas.o sort.c $(LIBS)

lu1fac.mexglx : lu1fac.c mi27lu.o mi15blas.o sort.c
	mex $(OPTIONS)  -L/common/gcc-4.2.2/bin/../lib/gcc/i686-pc-linux-gnu/4.2.2 -lgfortranbegin -lgfortran -lm -lgcc_s lu1fac.c -I$(INCLUDE) mi27lu.o mi15blas.o sort.c $(LIBS)

mi27lu.o : mi27lu.f
	g77 $(F77OPTS) -o mi27lu.o mi27lu.f

mi15blas.o : mi15blas.f
	g77 $(F77OPTS) -o mi15blas.o mi15blas.f

ZIP : 
	\rm -f lusolMex.zip
	zip    lusolMex.zip \
	   README makefile mi15blas.f mi27lu.f \
	   nrutil.h lu1fac.c sort.c \
	   lusolFactor.m luSOL.m lusolSet.m lusolTest.m \
	   illc1033.mat
	ls -lFh    --color=auto --time-style=locale *.zip

