#
# Unix or Linux script for testing LUSOL.
# Michael Saunders, Systems Optimization Laboratory,
# Stanford University.  saunders@stanford.edu
#
# 10 Sep 2002: First version for LUSOL distribution.
# 09 Oct 2003: lusoltest1.f written.
#              lusoltest2.f is not yet ready.


# System dependent defs go here.
# Main requirement: Select a compiler from one of the USE... choices below.

  OS        = $(shell uname -s)
  ARCH      = $(OS)-$(shell uname -m)
  MAKE      = make

# Set exactly one of these to yes

  USEabsoft = no
  USEg77    = yes
  USEifc    = no
  USElf95   = no
  USEsgi    = no
  USEsunos  = no
  USEgeneric= no

ifeq ($(USEabsoft),yes)    # LINUX with Absoft compiler
  FC        = f95
  FFLAGSdbg = -g
  FFLAGSopt = -O
  SHARED    = -X-shared
  LIBDIR    =
  LIBS      = -L$(ABSOFT) -lU77
endif

ifeq ($(USEg77),yes)       # LINUX with g77
  FC        = g77
  FFLAGSdbg = -g
  FFLAGSopt = -O
  SHARED    = -shared
  LIBDIR    =
  LIBS      = 
endif

ifeq ($(USEifc),yes)       # LINUX with Intel Fortran Compiler
  FC        = ifc
  SHARED    = --shared
  MFLAGS    = -v -f $(HOME)/matlab/mexopts90.sh
  FFLAGSdbg = -g
  FFLAGSopt = -O
# LIBDIR    = -L/usr/local/lib -lblas
  LIBS      = 
endif

ifeq ($(USElf95),yes)      # LINUX with Lahey compiler
# Lahey F95.
# For MINOS, have to use --nchk even with debug mode
# to accommodate passing real*8 z(lhs) as integer hs(nb), etc.
  FC        = lf95
  SHARED    = --shared
  MFLAGS    = -v -f $(HOME)/matlab/mexopts90.sh
# FFLAGSdbg = -g  --ap  --chk  --pca --nsav --tpp --trace --f95 --info
  FFLAGSdbg = -g  --ap --nchk  --pca --nsav --tpp --trace
  FFLAGSopt = -O --nap --nchk --npca --nsav --tpp --trace
  LIBDIR    = -L/usr/local/lib -lblas
  LIBS      = 
endif

ifeq ($(USEsgi),yes)       # SGI IRIX64
  FC77      = f77
  FC90      = f90
  FC        = f90
  MFLAGS    = -fortran -v
  FFLAGSdbg = -g -64
  FFLAGSopt = -O -64
  MFCFLA    = '-64 -freeform'
  FIXED     = -fixedform
  SHARED    = -shared      # <-- guessing!!
  LIBDIR    =
  LIBS      = 
endif

ifeq ($(USEsunos),yes)     # SunOS
  FC        = g77
  FFLAGSdbg = -g
  FFLAGSopt = -O
  SHARED    = -shared
  LIBDIR    =
  LIBS      = 
endif

ifeq ($(USEgeneric),yes)   # Generic Linux or Unix Fortran
  FC        = f90
  CC        = cc
  CFLAGS    = -O2
  MEX       = mex
  FFLAGSdbg = -g
  FFLAGSopt = -O
  SHARED    = -shared
  LIBDIR    =
  LIBS      = 
endif


# Set one of these
  FFLAGS  = $(FFLAGSdbg)
# FFLAGS  = $(FFLAGSopt)


# Clear suffix list, then define the ones we want
.SUFFIXES:
.SUFFIXES: .c .f .f90 .o

.f90.o:; ${FC} ${FFLAGS} -c -o $@ $<
.f.o:;   ${FC} ${FFLAGS} -c -o $@ $<
.c.o:;   $(CC) $(CFLAGS) -c -o $@ $<


LUSOL      = ../src/
filesLUSOL = ${LUSOL}lusol1.o \
             ${LUSOL}lusol2.o \
             ${LUSOL}lusol6a.o

files1     = ${filesLUSOL} mi15blas.o lusoltest1.o

files2     = ${filesLUSOL} mi15blas.o lusoltest2.o \
             gmscal.o hbsubs.o m1cpu.o dload.o dnrmx.o xerbla.o


lusoltest1:  ${files1}
	${FC} ${FFLAGS} -o $@ ${files1} ${LIBDIR}

lusoltest2:  ${files2}
	${FC} ${FFLAGS} -o $@ ${files2} ${LIBDIR}

clean:
	\rm -f *.o
