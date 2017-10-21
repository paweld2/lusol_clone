mxlu1fac.c and lu1fac.m together provides the cmex interface for
LU1FAC.  You need to compile the C program by turning on -DMATLAB and
-DYZHANG in the Makefile.  (-DMATLAB turns on matlab related memory
allocation/free routines; -DYZHANG turns on various fixes for various
bugs I found in the original LUSOL source code).

-- Yin Zhang
<yzhang@cs.utexas.edu>
UT-Austin
Dec 1, 2005
