   Notes on test programs in lusol/test
   Michael Saunders, SOL, Stanford University
   saunders@stanford.edu

10 Feb 2004:

lusoltest1.f is operational.
lusoltest2.f is not.

lusoltest1.f reads data from Afile.txt and bfile.txt
to obtain a sparse matrix A and a dense vector b.
It then uses lu1fac and lu6sol to factorize A and solve Ax = b
and judges whether the computed x gives a small residual r = b - Ax.

To test two sets of data obtained from Tao Gang (tao-gang@utulsa.edu),
proceed as follows.


FIRST EXAMPLE:

 make            # compiles lusoltest1.f and mi15blas.f
 cp -p A6805.txt Afile.txt
 cp -p b6805.txt bfile.txt
 ./lusoltest1

This should produce the following output:

 ==========
 lusoltest1
 ==========
 A  read successfully
 m      =    6805
 n      =    6805
 nnzero =   52660

 b  read successfully

 m        6805 =n        6805  Elems    52660  Amax   1.0E+09  Density   0.11
 Merit   294.9  lenL    85360  L+U     198622  Cmpressns    3  Incres  277.18
 Utri        0  lenU   113262  Ltol   5.0E+00  Umax   2.7E+09  Ugrwth 2.6E+00
 Ltri        0  dense1    421  Lmax   5.0E+00
 bump     6805  dense2    296  DUmax  2.4E+09  DUmin  1.3E+01  condU  1.8E+08

 bnorm =   1.3E+09
 rnorm =   6.5E-04
 xnorm =   3.1E+03
 rnorm is small.  The test seems successful.


SECOND EXAMPLE:

 cp -p A10009.txt Afile.txt
 cp -p b10009.txt bfile.txt
 ./lusoltest1

This should produce the following output:

 ==========
 lusoltest1
 ==========
 A  read successfully
 m      =   10009
 n      =   10009
 nnzero =   44590

 b  read successfully

 m       10009 =n       10009  Elems    44590  Amax   7.1E+04  Density   0.04
 Singular(m=n)  rank     7844  n-rank    2165  nsing     2165
 Merit   589.1  lenL   124779  L+U     189336  Cmpressns    3  Incres  324.62
 Utri        0  lenU    64557  Ltol   5.0E+00  Umax   2.2E+06  Ugrwth 3.1E+01
 Ltri     2232  dense1      0  Lmax   5.0E+00
 bump     7777  dense2      0  DUmax  7.1E+05  DUmin  1.1E-02  condU  6.5E+07

 bnorm =   7.0E+04
 rnorm =   2.1E+00
 xnorm =   5.7E+01
 rnorm is not very small.
 The LU factors may not be good enough.
 Try a smaller factol and/or Rook Pivoting

Note that A is determined to be very rank-deficient.
We cannot expect to get a small residual for Ax = b.
I have not been able to determine what Tao Gang expects
LUSOL to do in this case.

TPP (Threshold Partial Pivoting) with Ltol1 = 5.0 seems to determine
the rank correctly.
TRP (Threshold Rook Pivoting) with the same Ltol1 determines the same rank.
TRP with Ltol1 = 2.0 also determines the same rank (taking noticeably longer):

 ==========
 lusoltest1
 ==========
 A  read successfully
 m      =   10009
 n      =   10009
 nnzero =   44590

 b  read successfully

 m       10009 =n       10009  Elems    44590  Amax   7.1E+04  Density   0.04
 Singular(m=n)  rank     7844  n-rank    2165  nsing     2165
 MerRP   948.0  lenL   178666  L+U     261996  Cmpressns    7  Incres  487.57
 Utri        0  lenU    83330  Ltol   2.0E+00  Umax   7.1E+04  Ugrwth 1.0E+00
 Ltri     2138  dense1      0  Lmax   2.0E+00  Akmax  0.0E+00  Agrwth 0.0E+00
 bump     7871  dense2      0  DUmax  7.1E+04  DUmin  1.9E-03  condU  3.8E+07

 bnorm =   7.0E+04
 rnorm =   9.7E-01
 xnorm =   1.6E+01
 rnorm is not very small.
 The LU factors may not be good enough.
 Try a smaller factol and/or Rook Pivoting
