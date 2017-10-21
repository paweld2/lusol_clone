LUSOL maintains LU factors of a square or rectangular sparse matrix A.

The software for LUSOL is provided by SOL, Stanford University
under the terms of the Common Public License (CPL):
http://oss.software.ibm.com/developerworks/opensource/license-cpl.html

The src files lusol*.f and distribution are maintained by Michael Saunders.

09 Oct 2003: First set of files organized for distribution.  Includes
             matlab/* mex files contributed by Mike O'Sullivan.
06 Mar 2004: csrc/* files contributed by Kjell Eikland.
01 Dec 2005: cmex/* files contributed by Yin Zhang.
31 Jan 2008: Last update to this file.


Please send comments to
             Michael Saunders, SOL, Stanford University
             saunders@stanford.edu  650-723-1875
-----------------------------------------------------------------------------


DIRECTORIES

lusol           Top directory.
lusol/src       The LUSOL Fortran 77 source code.
lusol/matlab    cmex files by Mike O'Sullivan and some *.m files.
lusol/csrc      f77 -> Pascal -> C translations by Kjell Eikland.
lusol/cmex      The csrc files and mxlu1fac.c by Yin Zhang.
lusol/test      f77 test programs.


FILES

 1. LUSOL-overview.txt gives a general description of LUSOL.

 2. src/lusol.txt discusses the main f77 source files.

 3. src/mi27lu.f and src/sn27lu.f are the same.
    They contain the Factor, Solve, and Update routines used
    in the large-scale optimization packages MINOS and SNOPT.
    They are made by concatenating the following files in src:

      27HEAD.f
      lusol1.f
      lusol2.f
      lusol6a.f
      lusol7a.f
      lusol8a.f

 4. test/lusoltest1.f reads test data A, b and solves Ax = b.
    This is the test program to try first (for Fortran users).
    Check that the luparm(*) and parmlu(*) parameters have
    sensible values in lusoltest1.f, as documented near the
    top of lusol1.f.

    test/Makefile is intended for Unix and Linux systems
    to compile lusoltest1.f and associated source files.
    In case you are running Windows or some other operating system,
    the files needed are

      src/lusol1.f
      src/lusol2.f
      src/lusol6a.f
      test/mi15blas.f
      test/lusoltest1.f  (the main program)

 5. test/lusoltest2.f is not ready yet.

 6. test/mi25bfac.f and test/mi26bfac.f illustrate how the
    routines in src/mi27lu.f are used in MINOS and SNOPT for
    the following purposes:

    o Finding LU factors of a specified B.

    o Updating the factors when one column of B is replaced.

    o Factorizing the transpose of a larger matrix (B S)
      in order to find a column permutation P such that
          (B1 S1) = (B S)P
      has a reasonably well-conditioned square basis matrix B1
      and a well-conditioned null-space matrix operator

         Z = P[-B1\inv S1].
              [    I     ]


DOCUMENTATION

 1. The original reference for LUSOL is here:

      P. E. Gill, W. Murray, M. A. Saunders and M. H. Wright,
      Maintaining LU factors of a general sparse matrix,
      Linear Algebra and its Applications 88/89, 239-270 (1987).

 2. LUSOL and its use within MINOS and SNOPT are documented in
    sections 4 and 5 of the following paper:

      P. E. Gill, W. Murray and M. A. Saunders,
      SNOPT: An SQP algorithm for large-scale constrained optimization,
      SIGEST article, SIAM Review 47(1), 99-131 (2005).
      http://www.stanford.edu/group/SOL/papers/SNOPT-SIGEST.pdf

 3. LUSOL-overview.txt and src/lusol.txt contain general information.

 4. f77 comments in the *.f files give in-line documentation.
