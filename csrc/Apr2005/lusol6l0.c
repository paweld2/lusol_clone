
/* Create a row-based version of L0.
   This makes it possible to solve L0'x=h (btran) faster for sparse h,
   since we only run down the columns of L0' (rows of LO) for which
   the corresponding entry in h is non-zero. */
void LU1L0(LUSOLrec *LUSOL, LUSOLmat **mat, int *inform)
{
  int  K, L, LL, L1, L2, LEN, LENL0, NUML0, I, II;
  int  *lsumr;

  /* Check if there is anything worth doing and if so, initialize */
  if(mat == NULL)
    return;
  if(*mat != NULL)
    LUSOL_matfree(mat);
  NUML0 = LUSOL->luparm[LUSOL_IP_COLCOUNT_L0];
  LENL0 = LUSOL->luparm[LUSOL_IP_NONZEROS_L0];
  if((NUML0 == 0) || (LENL0 == 0))
    return;
  *mat = LUSOL_matcreate(LUSOL->n, LENL0);
  lsumr = (int *) malloc((LUSOL->n+1)*sizeof(*lsumr));
  if((*mat == NULL) || (lsumr == NULL)) {
    *inform = LUSOL_INFORM_NOMEMLEFT;
    return;
  }

  /* Compute non-zero counts by permuted row index */
  L2 = LUSOL->lena-LENL0;
  for(K = NUML0; K >= 1; K--) {
    LEN = LUSOL->lenc[K];
    L1 = L2+1;
    L2 += LEN;
    for(L = L1; L <= L2; L++) {
      I = LUSOL->indc[L];
      II = LUSOL->ipinv[I];
      (*mat)->vlen[II]++;
    }
  }

  /* Cumulate the row counts to get vector offsets; first row is leftmost  
     (Stick with Fortran array offset for consistency) */
  lsumr[1] = 1;
  for(K = 1; K < LUSOL->n; K++)
    lsumr[K+1] = lsumr[K] + (*mat)->vlen[K];

  /* Pack row lengths and update the row count parameter */
  L = 0;
  for(K = 1; K <= LUSOL->n; K++) {
    LEN = (*mat)->vlen[K];
    if(LEN > 0) {
      L++;
      (*mat)->vlen[L] = LEN;
    }
  }
  LUSOL->luparm[LUSOL_IP_ROWCOUNT_L0] = L;
  
  /* Map the matrix into row order by permuted index; first row is leftmost */
  L2 = LUSOL->lena-LENL0;
  for(K = 1; K <= NUML0; K++) {
    LEN = LUSOL->lenc[K];
    L1 = L2+1;
    L2 += LEN;
    for(L = L1; L <= L2; L++) {
      I = LUSOL->indc[L];
      II = LUSOL->ipinv[I];
      LL = lsumr[II]++;
      (*mat)->a[LL] = LUSOL->a[L];
      (*mat)->indr[LL] = LUSOL->indr[L];
      (*mat)->indc[LL] = I;
    }
  }

  /* Clean up */
  FREE(lsumr);
}

/* Solve L0' v = v based on row-based version of L0, constructed by LU1L0 */
void LU6L0T_v(LUSOLrec *LUSOL, LUSOLmat *mat, REAL V[], int NZidx[])
{
#ifdef DoTraceL0
  REAL TEMP;
#endif
  int  LEN, IPIV, K, L, L1, LENL0, NUML0;
  REAL SMALL;
  register REAL VPIV;
#ifdef LUSOLFastSolve
  REAL *aptr;
  int  *jptr;
#else
  int  J;
#endif

  NUML0 = LUSOL->luparm[LUSOL_IP_ROWCOUNT_L0];
  LENL0 = LUSOL->luparm[LUSOL_IP_NONZEROS_L0];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];

  /* Loop over the non-empty columns of L0' - from the end, going forward.
     Note that the diagonal is always 1, and not stored explicitly. */
  L1 = LENL0+1;
  for(K = NUML0; K > 0; K--) {
    LEN = mat->vlen[K];
    L = L1;
    L1 -= LEN;
    /* Get index and value of the corresponding entry of V[] */
    IPIV = mat->indc[L1];
    VPIV = V[IPIV];
    /* Only process the column of L0' if the value of V[] is non-zero */
    if(fabs(VPIV)>SMALL) {
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
      L--;
      for(aptr = mat->a+L, jptr = mat->indr+L;
          LEN > 0; LEN--, aptr--, jptr--)
        V[*jptr] += (*aptr) * VPIV;
#else
      for(; LEN > 0; LEN--) {
        L--;
        J = mat->indr[L];
#ifndef DoTraceL0
        V[J] += mat->a[L]*VPIV;
#else
        TEMP = V[J];
        V[J] += mat->a[L]*VPIV;
        printf("V[%3d] = V[%3d] + L[%d,%d]*V[%3d]\n", J, J, IPIV,J, IPIV);
        printf("%6g = %6g + %6g*%6g\n", V[J], TEMP, mat->a[L], VPIV);
#endif
      }
#endif
    }
  }

}

void print_L0T(LUSOLrec *LUSOL)
{
  int  I, J, K, L, L1, LEN, LENL0, NUML0, IPIV;
  LUSOLmat *mat = LUSOL->L0;
  REAL *denseL0 = (REAL*) calloc(LUSOL->m+1, (LUSOL->n+1)*sizeof(*denseL0));

  if(denseL0 == NULL)
    return;

  NUML0 = LUSOL->luparm[LUSOL_IP_COLCOUNT_L0];
  LENL0 = LUSOL->luparm[LUSOL_IP_NONZEROS_L0];

  L1 = 1;
  for(K = 1; K <= NUML0; K++) {
    LEN = mat->vlen[K];
    L = L1;
    L1 += LEN;
    IPIV = mat->indc[L];
    IPIV = LUSOL->ipinv[IPIV]; /* Undo row mapping */
    for(; LEN > 0; LEN--, L++) {
      J = mat->indr[L];
      denseL0[(LUSOL->n+1)*(J-1) + IPIV] = mat->a[L];
    }
  }

  for(I = 1; I <= LUSOL->n; I++) {
    for(J = 1; J <= LUSOL->m; J++)
      fprintf(stdout, "%10g", denseL0[(LUSOL->n+1)*(J-1) + I]);
    fprintf(stdout, "\n");
  }
  FREE(denseL0);
}
