
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   File  lusol6a
      lu6sol   lu6L     lu6Lt     lu6U     Lu6Ut   lu6LD   lu6chk
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   26 Apr 2002: lu6 routines put into a separate file.
   15 Dec 2002: lu6sol modularized via lu6L, lu6Lt, lu6U, lu6Ut.
                lu6LD implemented to allow solves with LDL' or L|D|L'.
   15 Dec 2002: Current version of lusol6a.f.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* ==================================================================
   lu6chk  looks at the LU factorization  A = L*U.
   If mode = 1, lu6chk is being called by lu1fac.
   (Other modes not yet implemented.)
   ------------------------------------------------------------------
   The important input parameters are

                  lprint = luparm(2)
                  keepLU = luparm(8)
                  Utol1  = parmlu(4)
                  Utol2  = parmlu(5)

   and the significant output parameters are

                  inform = luparm(10)
                  nsing  = luparm(11)
                  jsing  = luparm(12)
                  jumin  = luparm(19)
                  Lmax   = parmlu(11)
                  Umax   = parmlu(12)
                  DUmax  = parmlu(13)
                  DUmin  = parmlu(14)
                  and      w(*).

   Lmax  and Umax  return the largest elements in L and U.
   DUmax and DUmin return the largest and smallest diagonals of U
                   (excluding diagonals that are exactly zero).
   In general, w(j) is set to the maximum absolute element in
   the j-th column of U.  However, if the corresponding diagonal
   of U is small in absolute terms or relative to w(j)
   (as judged by the parameters Utol1, Utol2 respectively),
   then w(j) is changed to - w(j).
   Thus, if w(j) is not positive, the j-th column of A
   appears to be dependent on the other columns of A.
   The number of such columns, and the position of the last one,
   are returned as nsing and jsing.
   Note that nrank is assumed to be set already, and is not altered.
   Typically, nsing will satisfy      nrank + nsing = n,  but if
   Utol1 and Utol2 are rather large,  nsing > n - nrank   may occur.
   If keepLU = 0,
   Lmax  and Umax  are already set by lu1fac.
   The diagonals of U are in the top of A.
   Only Utol1 is used in the singularity test to set w(*).
   inform = 0  if  A  appears to have full column rank  (nsing = 0).
   inform = 1  otherwise  (nsing .gt. 0).
   ------------------------------------------------------------------
   00 Jul 1987: Early version.
   09 May 1988: f77 version.
   11 Mar 2001: Allow for keepLU = 0.
   17 Nov 2001: Briefer output for singular factors.
   05 May 2002: Comma needed in format 1100 (via Kenneth Holmstrom).
   06 May 2002: With keepLU = 0, diags of U are in natural order.
                They were not being extracted correctly.
   15 Dec 2002: Current version.
   ================================================================== */
void LU6CHK(LUSOLrec *LUSOL, int MODE, int LENA2, int *INFORM)
{
  MYBOOL KEEPLU;
  int    I, J, JSING, JUMIN, K, L, L1, L2, LENL, LPRINT, NDEFIC, NRANK, NSING;
  REAL   AIJ, DIAG, DUMAX, DUMIN, LMAX, UMAX, UTOL1, UTOL2;

  LPRINT = LUSOL->luparm[LUSOL_IP_PRINTLEVEL];
  KEEPLU = (MYBOOL) (LUSOL->luparm[LUSOL_IP_KEEPLU]!=0);
  NRANK = LUSOL->luparm[LUSOL_IP_RANK_U];
  LENL  = LUSOL->luparm[LUSOL_IP_NONZEROS_L];
  UTOL1 = LUSOL->parmlu[LUSOL_RP_SMALLDIAG_U];
  UTOL2 = LUSOL->parmlu[LUSOL_RP_EPSDIAG_U];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  LMAX  = ZERO;
  UMAX  = ZERO;
  NSING = 0;
  JSING = 0;
  JUMIN = 0;
  DUMAX = ZERO;
  DUMIN = LUSOL_BIGNUM;

#ifdef LUSOLFastClear
  MEMCLEAR(LUSOL->w, LUSOL->n + 1);
#else
  for(I = 1; I <= LUSOL->n; I++)
    LUSOL->w[I] = ZERO;
#endif

  if(KEEPLU) {
/*     --------------------------------------------------------------
        Find  Lmax.
       -------------------------------------------------------------- */
    for(L = (LENA2+1)-LENL; L <= LENA2; L++) {
      LMAX = MAX(LMAX,fabs(LUSOL->a[L]));
     }
/*     --------------------------------------------------------------
        Find Umax and set w(j) = maximum element in j-th column of U.
       -------------------------------------------------------------- */
    for(K = 1; K <= NRANK; K++) {
      I = LUSOL->ip[K];
      L1 = LUSOL->locr[I];
      L2 = (L1+LUSOL->lenr[I])-1;
      for(L = L1; L <= L2; L++) {
        J = LUSOL->indr[L];
        AIJ = fabs(LUSOL->a[L]);
        LUSOL->w[J] = MAX(LUSOL->w[J],AIJ);
        UMAX = MAX(UMAX,AIJ);
      }
    }
/*     --------------------------------------------------------------
        Negate w(j) if the corresponding diagonal of U is
        too small in absolute terms or relative to the other elements
        in the same column of  U.
        Also find DUmax and DUmin, the extreme diagonals of U.
       -------------------------------------------------------------- */
    for(K = 1; K <= LUSOL->n; K++) {
      J = LUSOL->iq[K];
      if(K>NRANK)
        DIAG = ZERO;
      else {
        I = LUSOL->ip[K];
        L1 = LUSOL->locr[I];
        DIAG = fabs(LUSOL->a[L1]);
        DUMAX = MAX(DUMAX,DIAG);
        if(DUMIN>DIAG) {
          DUMIN = DIAG;
          JUMIN = J;
        }
      }
      if(DIAG<=UTOL1 || DIAG<=UTOL2*LUSOL->w[J]) {
        NSING++;
        JSING = J;
        LUSOL->w[J] = -LUSOL->w[J];
      }
    }
    LUSOL->parmlu[LUSOL_RP_MAXMULT_L] = LMAX;
    LUSOL->parmlu[LUSOL_RP_MAXELEM_U] = UMAX;
  }
   else {
/*     --------------------------------------------------------------
        keepLU = 0.
        Only diag(U) is stored.  Set w(*) accordingly.
       -------------------------------------------------------------- */
    for(K = 1; K <= LUSOL->n; K++) {
      J = LUSOL->iq[K];
      if(K>NRANK)
        DIAG = ZERO;
      else {
/* !             diag   = abs( diagU(k) ) ! 06 May 2002: Diags are in natural order */
        DIAG = fabs(LUSOL->diagU[J]);
        LUSOL->w[J] = DIAG;
        DUMAX = MAX(DUMAX,DIAG);
        if(DUMIN>DIAG) {
          DUMIN = DIAG;
          JUMIN = J;
        }
      }
      if(DIAG<=UTOL1) {
        NSING++;
        JSING = J;
        LUSOL->w[J] = -LUSOL->w[J];
      }
    }
  }
/*     -----------------------------------------------------------------
        Set output parameters.
       ----------------------------------------------------------------- */
  if(JUMIN==0)
    DUMIN = ZERO;
  LUSOL->luparm[LUSOL_IP_SINGULARITIES]  = NSING;
  LUSOL->luparm[LUSOL_IP_SINGULARINDEX]  = JSING;
  LUSOL->luparm[LUSOL_IP_COLINDEX_DUMIN] = JUMIN;
  LUSOL->parmlu[LUSOL_RP_MAXELEM_DIAGU]  = DUMAX;
  LUSOL->parmlu[LUSOL_RP_MINELEM_DIAGU]  = DUMIN;
/*      The matrix has been judged singular. */
  if(NSING>0) {
    *INFORM = LUSOL_INFORM_LUSINGULAR;
    NDEFIC = LUSOL->n-NRANK;
    if(LPRINT>=LUSOL_MSG_SINGULARITY) {
      LUSOL_report(LUSOL, 0, "Singular(m%cn)  rank:%9d  n-rank:%8d  nsing:%9d\n",
                             relationChar(LUSOL->m, LUSOL->n),NRANK,NDEFIC,NSING);
    }
  }
/*      Exit. */
  LUSOL->luparm[LUSOL_IP_INFORM] = *INFORM;
}

/* ------------------------------------------------------------------
   lu6L   solves   L v = v(input).
   ------------------------------------------------------------------
   15 Dec 2002: First version derived from lu6sol.
   15 Dec 2002: Current version.
   ------------------------------------------------------------------ */
void LU6L(LUSOLrec *LUSOL, int *INFORM, REAL V[])
{
  int  IPIV, K, L, L1, LDUMMY, LEN, LENL, LENL0, NUML;
  REAL SMALL, VPIV, NUML0;
#ifdef LUSOLFastSolve
  REAL *aptr;
  int  *iptr, *jptr;
#else
  int  I, J;
#endif

  NUML0 = LUSOL->luparm[LUSOL_IP_COLCOUNT_L0];
  LENL0 = LUSOL->luparm[LUSOL_IP_NONZEROS_L0];
  LENL  = LUSOL->luparm[LUSOL_IP_NONZEROS_L];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  L1 = LUSOL->lena+1;
  for(K = 1; K <= NUML0; K++) {
    LEN = LUSOL->lenc[K];
    L = L1;
    L1 -= LEN;
    IPIV = LUSOL->indr[L1];
    VPIV = V[IPIV];
    if(fabs(VPIV)>SMALL) {
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
      L--;
      for(LDUMMY = 1, aptr = LUSOL->a+L, jptr = LUSOL->indc+L;
          LDUMMY <= LEN; LDUMMY++, aptr--, jptr--)
        V[*jptr] += (*aptr) * VPIV;
#else
      for(LDUMMY = 1; LDUMMY <= LEN; LDUMMY++) {
        L--;
        J = LUSOL->indc[L];
        V[J] += LUSOL->a[L]*VPIV;
      }
#endif
    }
  }
  L = (LUSOL->lena-LENL0)+1;
  NUML = LENL-LENL0;
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
  L--;
  for(LDUMMY = 1, aptr = LUSOL->a+L, iptr = LUSOL->indr+L, jptr = LUSOL->indc+L;
      LDUMMY <= NUML; LDUMMY++, aptr--, iptr--, jptr--) {
    if(fabs(V[*iptr])>SMALL)
      V[*jptr] += (*aptr) * V[*iptr];
  }
#else
  for(LDUMMY = 1; LDUMMY <= NUML; LDUMMY++) {
    L--;
    I = LUSOL->indr[L];
    if(fabs(V[I])>SMALL) {
      J = LUSOL->indc[L];
      V[J] += LUSOL->a[L]*V[I];
    }
  }
#endif
/*      Exit. */
  LUSOL->luparm[LUSOL_IP_INFORM] = *INFORM;
}

/* ==================================================================
   lu6LD  assumes lu1fac has computed factors A = LU of a
   symmetric definite or quasi-definite matrix A,
   using Threshold Symmetric Pivoting (TSP),   luparm(6) = 3,
   or    Threshold Diagonal  Pivoting (TDP),   luparm(6) = 4.
   It also assumes that no updates have been performed.
   In such cases,  U = D L', where D = diag(U).
   lu6LDL returns v as follows:

   mode
    1    v  solves   L D v = v(input).
    2    v  solves   L|D|v = v(input).
   ------------------------------------------------------------------
   15 Dec 2002: First version of lu6LD.
   15 Dec 2002: Current version.
   ================================================================== */
void LU6LD(LUSOLrec *LUSOL, int *INFORM, int MODE, REAL V[])
{
  int  IPIV, K, L, L1, LDUMMY, LEN, NUML0;
  REAL DIAG, SMALL, VPIV;
#ifdef LUSOLFastSolve
  REAL *aptr;
  int  *jptr;
#else
  int  J;
#endif

/*      Solve L D v(new) = v  or  L|D|v(new) = v, depending on mode.
        The code for L is the same as in lu6L,
        but when a nonzero entry of v arises, we divide by
        the corresponding entry of D or |D|. */
  NUML0 = LUSOL->luparm[LUSOL_IP_COLCOUNT_L0];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  L1 = LUSOL->lena+1;
  for(K = 1; K <= NUML0; K++) {
    LEN = LUSOL->lenc[K];
    L = L1;
    L1 -= LEN;
    IPIV = LUSOL->indr[L1];
    VPIV = V[IPIV];
    if(fabs(VPIV)>SMALL) {
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
      L--;
      for(LDUMMY = 1, aptr = LUSOL->a+L, jptr = LUSOL->indc+L;
          LDUMMY <= LEN; LDUMMY++, aptr--, jptr--)
        V[*jptr] += (*aptr) * VPIV;
#else
      for(LDUMMY = 1; LDUMMY <= LEN; LDUMMY++) {
        L--;
        J = LUSOL->indc[L];
        V[J] += LUSOL->a[L]*VPIV;
      }
#endif
/*      Find diag = U(ipiv,ipiv) and divide by diag or |diag|. */
      L = LUSOL->locr[IPIV];
      DIAG = LUSOL->a[L];
      if(MODE==2)
        DIAG = fabs(DIAG);
      V[IPIV] = VPIV/DIAG;
    }
  }
}

/* ==================================================================
   lu6Lt  solves   L'v = v(input).
   ------------------------------------------------------------------
   15 Dec 2002: First version derived from lu6sol.
   15 Dec 2002: Current version.
   ================================================================== */
void LU6LT(LUSOLrec *LUSOL, int *INFORM, REAL V[])
{
  int  IPIV, K, L, L1, L2, LEN, LENL, LENL0, NUML0;
  REAL SMALL, SUM;
#ifdef LUSOLFastSolve
  REAL *aptr;
  int  *iptr, *jptr;
#else
  int  I, J;
#endif

  NUML0 = LUSOL->luparm[LUSOL_IP_COLCOUNT_L0];
  LENL0 = LUSOL->luparm[LUSOL_IP_NONZEROS_L0];
  LENL  = LUSOL->luparm[LUSOL_IP_NONZEROS_L];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  L1 = (LUSOL->lena-LENL)+1;
  L2 = LUSOL->lena-LENL0;
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
  for(L = L1, aptr = LUSOL->a+L1, iptr = LUSOL->indr+L1, jptr = LUSOL->indc+L1;
      L <= L2; L++, aptr++, iptr++, jptr++) {
    if(fabs(V[*jptr])>SMALL)
      V[*iptr] += (*aptr) * V[*jptr];
  }
#else
  for(L = L1; L <= L2; L++) {
    J = LUSOL->indc[L];
    if(fabs(V[J])>SMALL) {
      I = LUSOL->indr[L];
      V[I] += LUSOL->a[L]*V[J];
    }
  }
#endif
  for(K = NUML0; K >= 1; K--) {
    LEN = LUSOL->lenc[K];
    SUM = ZERO;
    L1 = L2+1;
    L2 += LEN;
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
    for(L = L1, aptr = LUSOL->a+L1, jptr = LUSOL->indc+L1;
        L <= L2; L++, aptr++, jptr++)
      SUM += (*aptr) * V[*jptr];
#else
    for(L = L1; L <= L2; L++) {
      J = LUSOL->indc[L];
      SUM += LUSOL->a[L]*V[J];
    }
#endif
    IPIV = LUSOL->indr[L1];
    V[IPIV] += SUM;
  }
/*      Exit. */
  LUSOL->luparm[LUSOL_IP_INFORM] = *INFORM;
}

/* ==================================================================
   lu6U   solves   U w = v.          v  is not altered.
   ------------------------------------------------------------------
   15 Dec 2002: First version derived from lu6sol.
   15 Dec 2002: Current version.
   ================================================================== */
void LU6U(LUSOLrec *LUSOL, int *INFORM, REAL V[], REAL W[])
{
  int  I, J, K, KLAST, L, L1, L2, L3, NRANK, NRANK1;
  REAL RESID, SMALL, T;
#ifdef LUSOLFastSolve
  REAL *aptr;
  int  *jptr;
#endif

  NRANK = LUSOL->luparm[LUSOL_IP_RANK_U];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  NRANK1 = NRANK+1;
  RESID = ZERO;
/*      Find the first nonzero in v(1:nrank), counting backwards. */
  for(KLAST = NRANK; KLAST >= 1; KLAST--) {
    I = LUSOL->ip[KLAST];
    if(fabs(V[I])>SMALL)
      break;
  }
  for(K = KLAST+1; K <= LUSOL->n; K++) {
    J = LUSOL->iq[K];
    W[J] = ZERO;
  }
/*      Do the back-substitution, using rows 1:klast of U. */
  for(K = KLAST; K >= 1; K--) {
    I = LUSOL->ip[K];
    T = V[I];
    L1 = LUSOL->locr[I];
    L2 = L1+1;
    L3 = (L1+LUSOL->lenr[I])-1;
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
    for(L = L2, aptr = LUSOL->a+L2, jptr = LUSOL->indr+L2;
        L <= L3; L++, aptr++, jptr++)
      T -= (*aptr) * W[*jptr];
#else
    for(L = L2; L <= L3; L++) {
      J = LUSOL->indr[L];
      T -= LUSOL->a[L]*W[J];
    }
#endif
    J = LUSOL->iq[K];
    if(fabs(T)<=SMALL)
      W[J] = ZERO;
    else
      W[J] = T/LUSOL->a[L1];
  }
/*      Compute residual for overdetermined systems. */
  for(K = NRANK1; K <= LUSOL->m; K++) {
    I = LUSOL->ip[K];
    RESID += fabs(V[I]);
  }
/*      Exit. */
  if(RESID>ZERO)
    *INFORM = LUSOL_INFORM_LUSINGULAR;
  LUSOL->luparm[LUSOL_IP_INFORM]     = *INFORM;
  LUSOL->parmlu[LUSOL_RP_RESIDUAL_U] = RESID;
}

/* ==================================================================
   lu6Ut  solves   U'v = w.          w  is destroyed.
   ------------------------------------------------------------------
   15 Dec 2002: First version derived from lu6sol.
   15 Dec 2002: Current version.
   ================================================================== */
void LU6UT(LUSOLrec *LUSOL, int *INFORM, REAL V[], REAL W[])
{
  int  I, J, K, L, L1, L2, NRANK, NRANK1;
  REAL RESID, SMALL, T;
#ifdef LUSOLFastSolve
  REAL *aptr;
  int  *jptr;
#endif

  NRANK = LUSOL->luparm[LUSOL_IP_RANK_U];
  SMALL = LUSOL->parmlu[LUSOL_RP_ZEROTOLERANCE];
  *INFORM = LUSOL_INFORM_LUSUCCESS;
  NRANK1 = NRANK+1;
  RESID = ZERO;
  for(K = NRANK1; K <= LUSOL->m; K++) {
    I = LUSOL->ip[K];
    V[I] = ZERO;
  }
/*      Do the forward-substitution, skipping columns of U(transpose)
        when the associated element of w(*) is negligible. */
  for(K = 1; K <= NRANK; K++) {
    I = LUSOL->ip[K];
    J = LUSOL->iq[K];
    T = W[J];
    if(fabs(T)<=SMALL) {
      V[I] = ZERO;
      continue;
    }
    L1 = LUSOL->locr[I];
    T /= LUSOL->a[L1];
    V[I] = T;
    L2 = (L1+LUSOL->lenr[I])-1;
    L1++;
/*     ***** This loop could be coded specially. */
#ifdef LUSOLFastSolve
    for(L = L1, aptr = LUSOL->a+L1, jptr = LUSOL->indr+L1;
        L <= L2; L++, aptr++, jptr++)
      W[*jptr] -= T * (*aptr);
#else
    for(L = L1; L <= L2; L++) {
      J = LUSOL->indr[L];
      W[J] -= T*LUSOL->a[L];
    }
#endif
  }
/*      Compute residual for overdetermined systems. */
  for(K = NRANK1; K <= LUSOL->n; K++) {
    J = LUSOL->iq[K];
    RESID += fabs(W[J]);
  }
/*      Exit. */
  if(RESID>ZERO)
    *INFORM = LUSOL_INFORM_LUSINGULAR;
  LUSOL->luparm[LUSOL_IP_INFORM]     = *INFORM;
  LUSOL->parmlu[LUSOL_RP_RESIDUAL_U] = RESID;
}

/* ==================================================================
   lu6sol  uses the factorization  A = L U  as follows:
   ------------------------------------------------------------------
   mode
    1    v  solves   L v = v(input).   w  is not touched.
    2    v  solves   L'v = v(input).   w  is not touched.
    3    w  solves   U w = v.          v  is not altered.
    4    v  solves   U'v = w.          w  is destroyed.
    5    w  solves   A w = v.          v  is altered as in 1.
    6    v  solves   A'v = w.          w  is destroyed.

   If mode = 3,4,5,6, v and w must not be the same arrays.
   If lu1fac has just been used to factorize a symmetric matrix A
   (which must be definite or quasi-definite), the factors A = L U
   may be regarded as A = LDL', where D = diag(U).  In such cases,

   mode
    7    v  solves   A v = L D L'v = v(input).   w  is not touched.
    8    v  solves       L |D| L'v = v(input).   w  is not touched.

   ip(*), iq(*)      hold row and column numbers in pivotal order.
   lenc(k)           is the length of the k-th column of initial L.
   lenr(i)           is the length of the i-th row of U.
   locc(*)           is not used.
   locr(i)           is the start  of the i-th row of U.

   U is assumed to be in upper-trapezoidal form (nrank by n).
   The first entry for each row is the diagonal element
   (according to the permutations  ip, iq).  It is stored at
   location locr(i) in a(*), indr(*).

   On exit, inform = 0 except as follows.
     if(mode = 3,4,5,6 and if U (and hence A) is singular,)
     inform = 1 if there is a nonzero residual in solving the system
     involving U.  parmlu(20) returns the norm of the residual.
   ------------------------------------------------------------------
     July 1987: Early version.
   09 May 1988: f77 version.
   27 Apr 2000: Abolished the dreaded "computed go to".
                But hard to change other "go to"s to "if then else".
   15 Dec 2002: lu6L, lu6Lt, lu6U, lu6Ut added to modularize lu6sol.
   ================================================================== */
void LU6SOL(LUSOLrec *LUSOL, int MODE, REAL V[], REAL W[], int *INFORM)
{
/*      Solve  L v(new) = v. */
  if(MODE==LUSOL_SOLVE_Lv_v) {
    LU6L(LUSOL, INFORM,V);
/*      Solve  L'v(new) = v. */
  }
  else if(MODE==LUSOL_SOLVE_Ltv_v) {
    LU6LT(LUSOL, INFORM,V);
/*      Solve  U w = v. */
  }
  else if(MODE==LUSOL_SOLVE_Uw_v) {
    LU6U(LUSOL, INFORM,V,W);
/*      Solve  U'v = w. */
  }
  else if(MODE==LUSOL_SOLVE_Utv_w) {
    LU6UT(LUSOL, INFORM,V,W);
/*      Solve  A w      = v */
  }
  else if(MODE==LUSOL_SOLVE_Aw_v) {
/*      via    L v(new) = v */
    LU6L(LUSOL, INFORM,V);
/*       and    U w = v(new). */
    LU6U(LUSOL, INFORM,V,W);
/*      Solve  A'v = w */
  }
  else if(MODE==LUSOL_SOLVE_Atv_w) {
/*      via    U'v = w */
    LU6UT(LUSOL, INFORM,V,W);
/*      and    L'v(new) = v. */
    LU6LT(LUSOL, INFORM,V);
  }
  else if(MODE==LUSOL_SOLVE_Av_v) {
/*      Solve  LDv(bar) = v */
    LU6LD(LUSOL, INFORM,1,V);
/*      and    L'v(new) = v(bar). */
    LU6LT(LUSOL, INFORM,V);
  }
  else if(MODE==LUSOL_SOLVE_LDLtv_v) {
/*      Solve  L|D|v(bar) = v */
    LU6LD(LUSOL, INFORM,2,V);
/*      and    L'v(new) = v(bar). */
    LU6LT(LUSOL, INFORM,V);
  }
}

