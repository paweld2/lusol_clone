/* ============================================================= */
/* === MATLAB/mxlu1fac mexFunction ============================= */
/* ============================================================= */

/* $Header: /people/cvs/cvsroot/lusol/cmex/mxlu1fac.c,v 1.1 2006/02/28 02:39:35 saunders Exp $ */

#include <math.h>
#include <time.h>
#include "commonlib.h"
#include "myblas.h"
#include "lusol.h"
#include "mex.h"

void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )
{
  /* Declare variables */
  int Am,An,Anz,Pm,Pn,*iA,*jA,*Ap,*Ai,*Lp,*Li,*Up,*Ui,*occupied;
  double *Ax,*Px,*Lx,*Ux,*Qx,*Aij,aij;
  int pivot,m,inform,i,j,k,pa,paend;
  int numL0,Lnz,Unz,pos,nrank,len,*Lstart,*Llen;
  double Amax;
  mxArray *L,*P,*Upt,*Q;

  /* Storage for LUSOL */
  LUSOLrec *LUSOL = NULL;

  /* Check for proper number of input and output arguments. */    
  if (nrhs < 1 || nrhs > 4 || nrhs > 4) {
    mexErrMsgTxt("Usage: [L,U,p,q] = mxlusol(A,pivottype,tol,memscalar).");
  }

  /* Check data type of input argument. */
  if (!(mxIsDouble(prhs[0]))) {
    mexErrMsgTxt("Input arguments must be of type double.");
  }
  if (mxGetNumberOfDimensions(prhs[0]) != 2) {
    mexErrMsgTxt("First input argument must be two dimensional\n");
  }
 
  /* Get the size and pointers to input data. */
  Am  = mxGetM(prhs[0]);
  An  = mxGetN(prhs[0]);
  if (mxIsChar(prhs[0]) || !mxIsSparse(prhs[0]) || mxIsComplex(prhs[0])) {
    mexErrMsgTxt("First argument must be a sparse real matrix.");
  }
  Ax  = mxGetPr(prhs[0]);
  Ai  = mxGetIr(prhs[0]);
  Ap  = mxGetJc(prhs[0]);
  Anz = mxGetNzmax(prhs[0]);

  pivot = LUSOL_PIVOT_TRP;
  if (nrhs >= 2) {
    Pm  = mxGetM(prhs[1]);
    Pn  = mxGetN(prhs[1]);
    if (!(Pm == 1 && Pn == 1) ||
        !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
      mexErrMsgTxt("Second argument must be a real scalar.");
    }
    Px  = mxGetPr(prhs[1]);
    m = (int) *Px;
    if(m >= 0 && m <= LUSOL_PIVOT_MAX) {
      pivot = m;
    }
  }

  LUSOL = LUSOL_create(stdout, 0, pivot, 0);
  LUSOL->luparm[LUSOL_IP_SCALAR_NZA] = Anz < 500000 ? 10 : LUSOL_MULT_nz_a;
  LUSOL->luparm[LUSOL_IP_PRINTLEVEL] = LUSOL_MSG_NONE;
  
  if (nrhs >= 3) {
    Pm  = mxGetM(prhs[2]);
    Pn  = mxGetN(prhs[2]);
    if (!(Pm == 1 && Pn == 1) ||
        !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
      mexErrMsgTxt("Second argument must be a real scalar.");
    }
    Px  = mxGetPr(prhs[2]);
    Amax = *Px;
    if(Amax >= 1 && Amax <= 100) {
      LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij] = Amax;
    }
  }

  if (nrhs >= 4) {
    Pm  = mxGetM(prhs[3]);
    Pn  = mxGetN(prhs[3]);
    if (!(Pm == 1 && Pn == 1) ||
        !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
      mexErrMsgTxt("Second argument must be a real scalar.");
    }
    Px  = mxGetPr(prhs[3]);
    m = (int) *Px;
    if(m >= LUSOL_MULT_nz_a && m <= 100) {
      LUSOL->luparm[LUSOL_IP_SCALAR_NZA] = m;
    }
  }

  /* Create the arrays. */
  Aij = (REAL *) mxCalloc(Anz + BASE, sizeof(REAL));
  iA = (int *)   mxCalloc(Anz + BASE, sizeof(int));
  jA = (int *)   mxCalloc(Anz + BASE, sizeof(int));

  k = BASE;
  for (j = 0; j < An; j++) {
    paend = Ap[j+1];
    for (pa = Ap[j]; pa < paend; pa++) {
      i = Ai[pa];
      aij = Ax[pa];
      iA [k] = i+BASE;
      jA [k] = j+BASE;
      Aij[k] = aij;
      k++;        
    }
  }

  if (!LUSOL_assign(LUSOL, iA, jA, Aij, Anz, TRUE)) {
    mexErrMsgTxt("LUSOL failed due to insufficient memory.");
  }

  /* Free the space */
  mxFree(Aij);
  mxFree(iA);
  mxFree(jA);

  /* Factor  A = L U. */
  LU1FAC( LUSOL, &inform );
  if (inform > LUSOL_INFORM_SERIOUS) {
    mexErrMsgTxt(LUSOL_informstr(LUSOL, inform));
  }

  /* Extract vectors P and Q */
  nrank = LUSOL->luparm[LUSOL_IP_RANK_U];
  occupied = (int*) mxCalloc(Am+1, sizeof(int));
  P   = mxCreateDoubleMatrix(1,Am,0);
  Px  = mxGetPr(P);
  for (k = 0, pos = 0; k < Am; k++) {
    if (k < nrank) {
      j = LUSOL->ip[k+BASE]-BASE;
    } else {
      while(occupied[pos]) {
        pos++;
      }
      j = pos;
    }
    occupied[j] = TRUE;
    Px[k] = j+1;
  }
  mxFree(occupied);

  occupied = (int*) mxCalloc(An+1, sizeof(int));
  Q   = mxCreateDoubleMatrix(1,An,0);
  Qx  = mxGetPr(Q);
  for (k = 0, pos = 0; k < An; k++) {
    if (k < nrank) {
      j = LUSOL->iq[k+BASE]-BASE;
    } else {
      while(occupied[pos]) {
        pos++;
      }
      j = pos;
    }
    occupied[j] = TRUE;
    Qx[k] = j+1;
  }
  mxFree(occupied);

  /* Extract L */
  numL0 = LUSOL->luparm[LUSOL_IP_COLCOUNT_L0];
  
  Lnz = LUSOL->luparm[LUSOL_IP_NONZEROS_L];
  L   = mxCreateSparse(Am,Am,Lnz+Am,0);
  Lx  = mxGetPr(L);
  Li  = mxGetIr(L);
  Lp  = mxGetJc(L);

  Lstart = (int*) mxCalloc(Am+1, sizeof(int));
  Llen   = (int*) mxCalloc(Am+1, sizeof(int));
  for (i = 0, pos = LUSOL->lena; i < numL0; i++, pos-=len) {
    j         = LUSOL->indr[pos]-BASE;
    len       = LUSOL->lenc[i+1];
    Lstart[j] = pos;
    Llen[j]   = len;
  }

  for (j = 0, k = 0; j < Am; j++) {
    Lp[j] = k;
    for (i = 0, pos = Lstart[j]; i < Llen[j]; i++, k++, pos--) {
      Li[k] = LUSOL->indc[pos]-BASE;
      Lx[k] = LUSOL->a[pos];
    }
  }
  Lp[Am] = k;
  mxFree(Llen);
  mxFree(Lstart);

  /* Extract U(p,:)^T (i.e. transpose of U) */
  Unz = LUSOL->luparm[LUSOL_IP_NONZEROS_U];
  Upt = mxCreateSparse(An,Am,Unz,0);
  Ux  = mxGetPr(Upt);
  Ui  = mxGetIr(Upt);
  Up  = mxGetJc(Upt);
  for (k = 0, pos = BASE; k < Unz; k++, pos++) {
    Ui[k] = LUSOL->indr[pos]-BASE;
    Ux[k] = LUSOL->a[pos];
  }
  for (j = 0, k = 0; j < nrank; j++) {
    Up[j] = k;
    k += LUSOL->lenr[LUSOL->ip[j+BASE]];
  }
  for (j = nrank; j <= Am; j++) {
    Up[j] = k;
  }

  /* Return the results to MATLAB */
  if (nlhs >= 1) plhs[0] = L;
  if (nlhs >= 2) plhs[1] = Upt;
  if (nlhs >= 3) plhs[2] = P;
  if (nlhs >= 4) plhs[3] = Q;
}

