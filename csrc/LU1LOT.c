
Dear Michael,

I put together some code to create a row-sorted version of L0 following a
factorization (function preliminarily named LU1L0T).  I also created a
function (LU6L0T_v) to solve using the row-sorted L0.  The code compiles and
runs without exceptions being raised, but the results are invalid.  Since
this is going into the abysses of LUSOL logic / coding my guess is that you
can quickly spot what is wrong.  For your information, the "LUSOLmat"
structure is a record containing vectors for the NZ-values ("a"), the row
counts ("len") and the column indeces of the NZ values ("indx").  I envision
that all transposed versions of LU components will use this format.


Best regards
Kjell

/* Create a row-based version of L0 */
void LU1L0T(LUSOLrec *LUSOL, LUSOLmat **mat)
{
  int  K, L, L1, LEN, LENL0, NUML0, I, J;
  int  *lsumr;

  NUML0 = LUSOL->luparm[LUSOL_IP_COLCOUNT_L0];
  LENL0 = LUSOL->luparm[LUSOL_IP_NONZEROS_L0];
  if(*mat != NULL)
    LUSOL_matfree(mat);
  if(*mat == NULL)
    *mat = LUSOL_matcreate(LUSOL->m, LENL0);
  lsumr = (int *) malloc((LUSOL->m+1)*sizeof(*lsumr));
  
  /* Compute row non-zero counts */
  MEMCLEAR((*mat)->vlen, LUSOL->m+1);
  L1 = LUSOL->lena+1;
  for(K = 1; K <= NUML0; K++) {
    LEN = LUSOL->lenc[K];
    L = L1;
    L1 -= LEN;
    for(; LEN > 0; LEN--) {
      L--;
      J = LUSOL->indc[L];
      (*mat)->vlen[J]++;
    }
  }
  
  /* Cumulate the row counts to get vector offsets */
  lsumr[1] = 0;
  for(K = 1; K < LUSOL->m; K++)
    lsumr[K+1] = lsumr[K] + (*mat)->vlen[K];
  
  /* Map the matrix into row order */
  MEMCLEAR((*mat)->vlen, LUSOL->m+1);
  L1 = LUSOL->lena+1;
  for(K = 1; K <= NUML0; K++) {
    LEN = LUSOL->lenc[K];
    L = L1;
    L1 -= LEN;
    I = LUSOL->indr[L1];
    for(; LEN > 0; LEN--) {
      L--;
      J = LUSOL->indc[L];
      (*mat)->vlen[J]++;
      J = lsumr[J]+(*mat)->vlen[J];
      (*mat)->a[J] = LUSOL->a[L];
      (*mat)->indx[J] = I;
    }
  }
  
  /* Clean up */
  FREE(lsumr);
}
/* Solve L0T v = v based on a row-based version of L0 */
void LU6L0T_v(LUSOLrec *LUSOL, LUSOLmat *mat, REAL V[], int NZidx[])
{
  int  LEN, IPIV, K, L, L1, LENL0, NUML0;
  REAL SMALL;
  register REAL VPIV;
#ifdef LUSOLFastSolve
  REAL *aptr;
  int  *jptr;
#else
  int  I, J;
#endif

  NUML0 = LUSOL->m;
  LENL0 = LUSOL->luparm[LUSOL_IP_NONZEROS_L0];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  L1 = LENL0+1;
  for(K = 1; K <= NUML0; K++) {
    LEN = mat->vlen[K];
    if(LEN == 0)
      continue;
    L = L1;
    L1 -= LEN;
    IPIV = mat->indx[L1];
    VPIV = V[IPIV];
    if(fabs(VPIV)>SMALL) {
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
      L--;
      for(aptr = mat->a+L, jptr = mat->indx+L;
          LEN > 0; LEN--, aptr--, jptr--)
        V[*jptr] += (*aptr) * VPIV;
#else
      for(; LEN > 0; LEN--) {
        L--;
        J = mat->indx[L];
        V[J] += mat->a[L]*VPIV;
      }
#endif
    }
  }
}


