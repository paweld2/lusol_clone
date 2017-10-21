/* lu1fac.c : written by M. J. Kim-O'Sullivan 20.1.99
   This is a CMEX file that calls LU1FAC, a sparse LU factorization written in 
   FORTRAN by Mike Saunders (MS). The MATLAB call is

   [L, U, p, q, luparm, parmlu] = lu1fac(A, luparm, parmlu)

   where A is a sparse MATLAB matrix, luparm and parmlu are LU
   settings, L and U are the sparse factors, and p and q are the row
   and column permutation vectors, respectively. */

#define INTEGER int

#include <string.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
 
  /* Inputs for LU1FAC */
  INTEGER m, n, nelem, lena;
  INTEGER luparm[30];
  double parmlu[30], *a, *w;
  INTEGER *indc, *indr, *ip, *iq,
    *lenc, *lenr,
    *iploc, *iqloc, *ipinv, *iqinv;
  INTEGER *locc, *locr;
  INTEGER inform;
  
  /* Local variables */
  mxArray *mp;
  double *pr;
  int *cnt, col, i, *inz, *ir, j, *jc, k, k1, l, l1, len, nz, start, stop;

  if (nrhs == 3) {

    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    pr = mxGetPr(prhs[0]);
    ir = mxGetIr(prhs[0]);
    jc = mxGetJc(prhs[0]);

    nelem = *(jc + n);
    /* printf("The matrix has dimensions (%d, %d), and %d elements.\n",
       m, n, nelem); */

    /* lena > max(2 * nelem, 10 * m, 10 * n, 500000) */
    lena = 500000;
    if (lena < 2 * nelem) lena = 2 * nelem;
    if (lena < 10 * m) lena = 10 * m;
    if (lena < 10 * n) lena = 10 * n;
    lena = lena + 1;
    /* printf("The length of a is %d.\n", lena); */
      
    /* Allocate all the LU1FAC arrays */
    a = (double *)calloc(lena, sizeof(double));
    indc = (INTEGER *)calloc(lena, sizeof(INTEGER));
    indr = (INTEGER *)calloc(lena, sizeof(INTEGER));
    ip = (INTEGER *)calloc(m, sizeof(INTEGER));
    iq = (INTEGER *)calloc(n, sizeof(INTEGER));
    
    lenc = (INTEGER *)calloc(n, sizeof(INTEGER));
    lenr = (INTEGER *)calloc(m, sizeof(INTEGER));
    locc = (INTEGER *)calloc(n, sizeof(INTEGER));
    locr = (INTEGER *)calloc(m, sizeof(INTEGER));
      
    iploc = (INTEGER *)calloc(n, sizeof(INTEGER));
    iqloc = (INTEGER *)calloc(m, sizeof(INTEGER));
    ipinv = (INTEGER *)calloc(m, sizeof(INTEGER));
    iqinv = (INTEGER *)calloc(n, sizeof(INTEGER));
    w = (double *)calloc(n, sizeof(double));

    /* Check if the allocations were successful */
    if (a && indc && indr && ip && iq && lenc && lenr && locc && locr &&
	iploc && iqloc && ipinv && iqinv && w) {
      
      /* Read the MATLAB matrix into the LU1FAC arrays */
      memcpy(a, (double *)pr, nelem * sizeof(double));
      nz = 0;
      for (i = 0; i < n; i++) {
	start = jc[i];
	stop  = jc[i + 1];
	if (start != stop)
	  for (j = start; j < stop; j++) {
	    indc[nz] = ir[j] + 1;
	    indr[nz] = i + 1;
	    nz++;
	  }
      }

      /* Extract luparm and parmlu */
      pr = mxGetPr(prhs[1]);
      for (i = 0; i < 30; i++)
	luparm[i] = (INTEGER)pr[i];
      
      pr = mxGetPr(prhs[2]);
      memcpy(parmlu, (double *)pr, 30 * sizeof(double));
      
      /* printf("Calling MINOS lu1fac.\n"); */
      /* Call LU1FAC */
      lu1fac_(&m, &n, &nelem, &lena, luparm, parmlu,
	      a, indc, indr, ip, iq,
	      lenc, lenr, locc, locr,
	      iploc, iqloc, ipinv, iqinv, w, &inform);
      
      if (inform <= 1) {
	if (luparm[7] > 0) {

	  /* Extract L from a, the diagonal of L is 1 and then the
	     rest is given by the non-trivial entries, for a total of
	     luparm[20] + m non-zeros */
	  mp = mxCreateSparse(m, m, luparm[20] + m, mxREAL);
	  pr = mxGetPr(mp);
	  ir = mxGetIr(mp);
	  jc = mxGetJc(mp);
	  cnt = (int *)calloc(n, sizeof(int));
	  if (!cnt) mexErrMsgTxt("Couldn't allocate memory");
	  jc[0] = 0;
	  for (k = 1; k <= n; k++) {
	    jc[k] = jc[k - 1] + 1;
	    l1 = lena - 1;
	    k1 = 0;
	    while ((k1 < luparm[19]) && (k != indr[l1]))
	      l1 = l1 - lenc[k1++];
	    if (k1 < luparm[19]) {
	      len = lenc[k1];
	      jc[k] += len;
	      for (l = 0; l < len; l++) {
		pr[jc[k - 1] + cnt[k - 1]] = -a[l1 - l];
		ir[jc[k - 1] + cnt[k - 1]++] = indc[l1 - l] - 1;
	      }
	      pr[jc[k - 1] + cnt[k - 1]] = 1.0;
	      ir[jc[k - 1] + cnt[k - 1]++] = k - 1;
	      /* Numerical recipes quicksort sorts arr[1..n], so pass the
		 address of the element before the start of the array
		 - MJKO. 02-03-99 */
	      sort(cnt[k - 1], &ir[jc[k - 1] - 1], &pr[jc[k - 1] - 1]);
	    } else {
	      pr[jc[k - 1] + cnt[k - 1]] = 1.0;
	      ir[jc[k - 1] + cnt[k - 1]++] = k - 1;
	    }
	  }
	  /* Replace the current L with the new one */
	  mxDestroyArray(plhs[0]);
	  plhs[0] = mp;
	  /* Free temporary pointers */
	  free((void *)cnt);
	  
	  /* Extract U from a */
	  mp = mxCreateSparse(m, n, luparm[21], mxREAL);
	  pr = mxGetPr(mp);
	  ir = mxGetIr(mp);
	  jc = mxGetJc(mp);
	  inz = (int *)calloc(n, sizeof(int));
	  if (!inz) mexErrMsgTxt("Couldn't allocate memory");
	  jc[0] = 0;
	  for (k = 0; k < luparm[15]; k++) {
	    i = ip[k];
	    inz[i - 1] = 1;
	    for (j = locr[i - 1]; j < (locr[i - 1] + lenr[i - 1]); j++)
	      jc[indr[j - 1]]++;
	  }
	  for (k = 0; k < n; k++) jc[k + 1] += jc[k];
	  
	  cnt = (int *)calloc(n, sizeof(int));
	  if (!cnt) mexErrMsgTxt("Couldn't allocate memory");
	  for (k = 0; k < m; k++)
	    if (inz[k]) {
	      for (j = locr[k]; j < (locr[k] + lenr[k]); j++) {
		col = indr[j - 1] - 1;
		pr[jc[col] + cnt[col]] = a[j - 1];
		ir[jc[col] + cnt[col]++] = k;
	      }
	    }
	  /* Replace the current U with the new one */
	  mxDestroyArray(plhs[1]);
	  plhs[1] = mp;
	  /* Free temporary pointers */
	  free((void *)inz);
	  free((void *)cnt);
	} /* End if keepLU */

	mp = mxCreateDoubleMatrix(1, m, mxREAL);    /* Extract p */
	pr = mxGetPr(mp);
	for (i = 0; i < m; i++)
	  pr[i] = ip[i];
	/* Replace the current matrix with the new one */
	mxDestroyArray(plhs[2]);
	plhs[2] = mp;
	
	mp = mxCreateDoubleMatrix(1, n, mxREAL);    /* Extract q */
	pr = mxGetPr(mp);
	for (i = 0; i < n; i++)
	  pr[i] = iq[i];
	/* Replace the current matrix with the new one */
	mxDestroyArray(plhs[3]);
	plhs[3] = mp;

	mp = mxCreateDoubleMatrix(1, 30, mxREAL);   /* Return luparm */
	pr = mxGetPr(mp);
	for (i = 0; i < 30; i++)
	  pr[i] = (double)luparm[i];
	/* Replace the current matrix with the new one */
	mxDestroyArray(plhs[4]);
	plhs[4] = mp;
      
	mp = mxCreateDoubleMatrix(1, 30, mxREAL);   /* Return parmlu */
	pr = mxGetPr(mp);
	memcpy(pr, (double *)parmlu, 30 * sizeof(double));
	/* Replace the current matrix with the new one */
	mxDestroyArray(plhs[5]);
	plhs[5] = mp;

      }
      else {
	printf("Inform returned from lu1fac = %d.\n", inform);
	if (inform == 7)
	  printf("Minimum value for lena should be %d.\n", luparm[12]);
	mexErrMsgTxt("LU1FAC returned a error");
      }
    }
    else mexErrMsgTxt("Couldn't allocate memory");
    
    /* Free all the LU1FAC variables */
    free((void *)a);
    free((void *)indc);
    free((void *)indr);
    free((void *)ip);
    free((void *)iq);
    free((void *)lenc);
    free((void *)lenr);
    free((void *)locc);
    free((void *)locr);
    free((void *)iploc);
    free((void *)iqloc);
    free((void *)ipinv);
    free((void *)iqinv);
    free((void *)w);
  }
  else mexErrMsgTxt("Invalid number of inputs");
  
} /* End of lu1fac mex-function */
