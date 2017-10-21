
/*     This program solves a sparse linear system Ax = b using LUSOL. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "commonlib.h"
#include "myblas.h"
#include "lusol.h"
#include "lusolio.h"
#include "lusolmain.h"

#ifdef YZHANG
/* XXX: Note that _strupr() is not available on many unix platforms */
void _strupr_(char *s)
{
  int  i;
  char c;
  for (i = 0; i < strlen(s); i++) {
    c = s[i];
    if (c <= 'z' && c >= 'a') {
      s[i] = c - 'a' + 'A';
    }
  }
}
#endif

MYBOOL getFileName(char *filename, char *test)
{
  MYBOOL status;
  status = (MYBOOL) (('-' != test[0]) && (strlen(test) > 1));
  if(status)
    strcpy(filename, test);
  return(status);
}
MYBOOL isNum(char val)
{
  int ord;
  ord = (int) val - 48;
  return( (MYBOOL) ((ord >= 0) && (ord <= 9)) );
}

int main( int argc, char *argv[], char *envp[] )
{
/* Output device */
  FILE *outunit = stdout;

/* Overall dimensions allocated */
  int    maxm = MAXROWS, maxn = MAXCOLS, maxnz = MAXNZ,
         replace = 0, randcol = 0;
  MYBOOL ftran = TRUE;

/* Storage for A, b */
#ifdef YZHANG
  /* XXX: initialize to shut gcc up */
  REAL   *Aij=0, *b=0, *xexact=0;
  int    *iA=0, *jA=0;
#else
  REAL   *Aij, *b, *xexact;
  int    *iA, *jA;
#endif
  
/* Storage for LUSOL */
  LUSOLrec *LUSOL = NULL;

/* Define local storage variables */
  int  i     , inform, j     , k     , i1   ,
       m     , useropt, lenb, lenx,
       n     , nelem , nnzero;
  REAL Amax  , test  ,
       bnorm , rnorm , xnorm,
#ifndef YZHANG    
       *rhs  , *r    , *x;
#else
       *rhs=0, *r=0  , *x=0;
#endif
  char fileext[50], filename[255], blasname[255];
  MYBOOL printsolution = FALSE, success = FALSE;

/* Check if we have input parameters */
  useropt = argc;
  if(useropt < 2) {
    printf("LUSOL v2.2 - Linear equation solver - Development edition.\n");
    printf("Usage:   lusol [-p <type>] [-t <tolerance>] [-m <memscalar>]\n");
    printf("               [-s] <MatrixFile> [<RhsFile>]\n");
    printf("Options: -p <type>        <0..3>    Selects pivot type.\n");
    printf("         -t <tolerance>   <1..100>  Selects diagonal tolerance.\n");
    printf("         -m <memscalar>   <%d..100> Memory allocation scalar.\n", LUSOL_MULT_nz_a);
    printf("         -b               Solves using btran (rather than ftran).\n");
    printf("         -r <times>       Randomly replace a column and resolve.\n");
    printf("         -s               Toggles display of the computed solution.\n");
    printf("         -blas <lib>      Activates exteral optimized BLAS library.\n");
    printf("Formats: Conventional RCV .TXT, MatrixMarket .MXT, or Harwell-Boeing .RUA\n");
    printf("Author:  Michael A. Saunders (original Fortran version)\n");
    printf("         Kjell Eikland       (modified C version)\n");
    return(0);
  }

/* Create the LUSOL object and set user options */
  LUSOL = LUSOL_create(outunit, 0, LUSOL_PIVOT_TPP, 0);
  LUSOL->luparm[LUSOL_IP_SCALAR_NZA] = 10;
  i = 1;
  n = 0;
  filename[0] = '\0';
  blasname[0] = '\0';
  while((n == 0) && (i < argc)) {
    if(strcmp("-p", argv[i]) == 0) {
      i1 = i+1;
#ifdef YZHANG
      /* XXX: bug fix */
      if((i1 < argc) && isNum(argv[i1][0])) {
#else
      if((i1 < argc) && isNum(argv[i1][1])) {
#endif
        i = i1;
        m = atoi(argv[i]);
        if(m < 0 || m > LUSOL_PIVOT_MAX)
          continue;
        LUSOL->luparm[LUSOL_IP_PIVOTTYPE] = m;
      }
    }
    else if(strcmp("-t", argv[i]) == 0) {
      i1 = i+1;
#ifdef YZHANG      
      /* XXX: bug fix */
      if((i1 < argc) && isNum(argv[i1][0])) {
#else
      if((i1 < argc) && isNum(argv[i1][1])) {
#endif
        i = i1;
        Amax = atof(argv[i]);
        if(Amax < 1 || Amax > 100)
          continue;
        LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij] = Amax;
      }
    }
    else if(strcmp("-m", argv[i]) == 0) {
      i1 = i+1;
#ifdef YZHANG      
      /* XXX: bug fix */
      if((i1 < argc) && isNum(argv[i1][0])) {
#else
      if((i1 < argc) && isNum(argv[i1][1])) {
#endif
        i = i1;
        m = atoi(argv[i]);
        if(m < LUSOL_MULT_nz_a || m > 100)
          continue;
        LUSOL->luparm[LUSOL_IP_SCALAR_NZA] = atoi(argv[i]);
      }
    }
    else if(strcmp("-s", argv[i]) == 0)
      printsolution = TRUE;
    else if(strcmp("-b", argv[i]) == 0)
      ftran = FALSE;
    else if(strcmp("-r", argv[i]) == 0) {
      i1 = i+1;
#ifdef YZHANG      
      /* XXX: bug fix */
      if((i1 < argc) && isNum(argv[i1][0])) {
#else
      if((i1 < argc) && isNum(argv[i1][1])) {
#endif
        i = i1;
        m = atoi(argv[i]);
        if(m < 0 || m > 10)
          continue;
      }
      else
        m = 1;
      srand((unsigned) time( NULL ));
      replace = 2*m;
    }
    else if(strcmp("-blas", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && getFileName(blasname, argv[i1])) {
        if(!load_BLAS(blasname))
          fprintf(outunit, "Could not load external BLAS library %s'\n", blasname);
        i = i1;
      }
      else
        fprintf(outunit, "Ignoring incomplete parameter %d '%s'\n", i, argv[i]);
    }
    else {
      if(getFileName(filename, argv[i])) {
        useropt = i;
        break;
      }
      else
        fprintf(outunit, "Ignoring unknown parameter %d '%s'\n", i, argv[i]);
    }
    i++;
  }

/* Obtain file extension and see if we must estimate matrix data size */
  strcpy(fileext, strchr(argv[useropt], '.'));
#ifdef YZHANG
  /* XXX: Note that _strupr() is not available on many unix platforms. */
  _strupr_(fileext);
#else
  _strupr(fileext);
#endif
  
  /* Read conventional text file format */
  if(strcmp(fileext, ".TXT") == 0) {
    if(!ctf_size_A(filename, &maxm, &maxn, &maxnz))
      goto x900;
  }
  /* Read MatrixMarket file format */
  else if(strcmp(fileext, ".MTX") == 0) {
    if(!mmf_size_A(filename, &maxm, &maxn, &maxnz))
      goto x900;
  }
  /* Read Harwell-Boeing file format */
  else if(strcmp(fileext, ".RUA") == 0) {
    if(!hbf_size_A(filename, &maxm, &maxn, &maxnz))
      goto x900;
  }
  else {
    fprintf(outunit, "Unrecognized matrix file extension %s\n", fileext);
    goto x900;
  }

/* Create the arrays */

  Aij = (REAL *) calloc(maxnz + BASE, sizeof(REAL));
  iA = (int *)   calloc(maxnz + BASE, sizeof(int));
  jA = (int *)   calloc(maxnz + BASE, sizeof(int));
  if(ftran)
    lenb = maxm;
  else
    lenb = maxn;
  b   = (REAL *) calloc(lenb+BASE, sizeof(REAL));
  rhs = (REAL *) calloc(lenb+BASE, sizeof(REAL));

  if(ftran)
    lenx = maxn;
  else
    lenx = maxm;
  xexact = (REAL *) calloc(lenx+BASE, sizeof(REAL));
  r = (REAL *) calloc(lenx+BASE, sizeof(REAL));
  x = (REAL *) calloc(lenx+BASE, sizeof(REAL));

/* -----------------------------------------------------------------
   Read data files.
   ----------------------------------------------------------------- */
  fprintf(stdout, "\n===================================\n");
  fprintf(stdout,   "LUSOL v2.2 - Linear equation solver");
  fprintf(stdout, "\n===================================\n");

/* -----------------------------------------------------------------
   Read data for A
   ----------------------------------------------------------------- */
  /* Read conventional text file format */
  if(strcmp(fileext, ".TXT") == 0) {
    if(!ctf_read_A(argv[useropt], maxm, maxn, maxnz,
                   &m, &n, &nnzero, iA, jA, Aij))
      goto x900;
  }
  /* Read MatrixMarket file format */
  else if(strcmp(fileext, ".MTX") == 0) {
    if(!mmf_read_A(argv[useropt], maxm, maxn, maxnz,
                   &m, &n, &nnzero, iA, jA, Aij))
      goto x900;
  }
  /* Read Harwell-Boeing file format */
  else if(strcmp(fileext, ".RUA") == 0) {
    if(!hbf_read_A(argv[useropt], maxm, maxn, maxnz,
                   &m, &n, &nnzero, iA, jA, Aij))
      goto x900;
  }
  else {
    fprintf(outunit, "Error: Unrecognized matrix file extension %s\n", fileext);
    goto x900;
  }

/* -----------------------------------------------------------------
   Read data for b
   ----------------------------------------------------------------- */
  /* Handle Harwell-Boeing case where the file contains a RHS */
  i = strcmp(fileext, ".RUA");

  if((useropt == argc) && (i != 0)) {
    i = (int) ((1.0*rand()/RAND_MAX)*(m-1));
    b[i+1] = 1.0;
  }
  else {
    if(i != 0)
      useropt++;
    strcpy(fileext, strchr(argv[useropt], '.'));
#ifdef YZHANG
    /* XXX: Note that _strupr() is not available on many unix platforms. */
    _strupr_(fileext);
#else
    _strupr(fileext);
#endif
    
    /* Read conventional text file format */
    if(strcmp(fileext, ".TXT") == 0) {
      if(!ctf_read_b(argv[useropt], lenb, b))
        goto x900;
    }
    /* Read MatrixMarket file format */
    else if(strcmp(fileext, ".MTX") == 0) {
      if(!mmf_read_b(argv[useropt], lenb, b))
        goto x900;
    }
  /* Read Harwell-Boeing file format */
    else if(strcmp(fileext, ".RUA") == 0) {
      if(!hbf_read_b(argv[useropt], lenb, b))
        goto x900;
    }
    else {
      fprintf(outunit, "Error: Unrecognized vector file extension %s\n", fileext);
      goto x900;
    }
  }
  success = TRUE;

/* -----------------------------------------------------------------
   Show data on input
   ----------------------------------------------------------------- */
  fprintf(outunit, "\nData read from:\n%s\n", filename);
  test = 100.0*nnzero/(m*n);
  fprintf(outunit, "Rows = %d   Columns = %d   Non-zeros = %d  Density =%5.1f%%\n",
                   m, n, nnzero, test);

/* -----------------------------------------------------------------
   Load A into (a, indc, indr).
   ----------------------------------------------------------------- */
#ifdef LegacyTesting
  LUSOL->luparm[LUSOL_IP_SCALAR_NZA] = LUSOL_MULT_nz_a;
/*  LUSOL->luparm[LUSOL_IP_PIVOTTYPE] = 1; */
  LUSOL_sizeto(LUSOL, MAXROWS, MAXCOLS, MAXNZ*LUSOL_MULT_nz_a);
#endif

  if(!LUSOL_assign(LUSOL, iA, jA, Aij, nnzero, TRUE)) {
    fprintf(outunit, "Error: LUSOL failed due to insufficient memory.\n");
    goto x900;
  }

/* ------------------------------------------------------------------
   Factor  A = L U.
   ------------------------------------------------------------------ */
  nelem = nnzero;
  LU1FAC( LUSOL, &inform );
  if (inform > LUSOL_INFORM_SERIOUS) {
    fprintf(outunit, "Error:\n%s\n", LUSOL_informstr(LUSOL, inform));
    goto x900;
  }
  /* Get the largest element in A; we use it below as an estimate
     of ||A||_inf, even though it isn't a proper norm. */
  Amax = LUSOL->parmlu[LUSOL_RP_MAXELEM_A];

/* ------------------------------------------------------------------
   SOLVE  A x = b.
   Save b first because lu6sol() overwrites rhs when computing x.
   ------------------------------------------------------------------ */
Resolve:
#if 1
  MEMCOPY(x, b, lenb+BASE);
  if(ftran)
    inform = LUSOL_ftran(LUSOL, x, NULL, FALSE);
  else
    inform = LUSOL_btran(LUSOL, x, NULL);
#else
  MEMCOPY(rhs, b, lenb+BASE);
  if(ftran)
    LU6SOL( LUSOL, LUSOL_SOLVE_Aw_v, rhs, x, NULL, &inform );
  else
    LU6SOL( LUSOL, LUSOL_SOLVE_Atv_w, x, rhs, NULL, &inform );
#endif
  if (inform > LUSOL_INFORM_SERIOUS) {
    fprintf(outunit, "Error:\n%s\n", LUSOL_informstr(LUSOL, inform));
    goto x900;
  }
  if(printsolution) {
    blockWriteREAL(outunit, "\nSolution vector", x, 1, lenb);
  }

/* ------------------------------------------------------------------
   Set r = b - Ax.
   Find norm of r and x.
   ------------------------------------------------------------------ */
  MEMCOPY(r, b, lenb+BASE);
  for (k = 1; k <= nnzero; k++) {
    i    = iA[k];
    j    = jA[k];
    if(ftran)
      r[i] -= Aij[k]*x[j];
    else
      r[j] -= Aij[k]*x[i];
  }
  bnorm  = dnormi( lenb, b );
  rnorm  = dnormi( lenb, r );
  xnorm  = dnormi( lenx, x );

/* ------------------------------------------------------------------
   Report the findings.
   ------------------------------------------------------------------ */
  if(randcol > 0)
    fprintf(outunit, "\n\nColumn %d was %s\n",
                      randcol,
                      (mod(replace,2) == 1 ? "replaced with random data" : "restored"));

  fprintf(outunit, "\nLU size statistics (%d reallocations):\n",
                   LUSOL->expanded_a);
  test = LUSOL->luparm[LUSOL_IP_NONZEROS_U0]+LUSOL->luparm[LUSOL_IP_NONZEROS_ROW];
  fprintf(outunit, "L0-size = %d   U0-size = %d   LU-nonzeros = %d   Fill-in = %.1fx\n",
                   LUSOL->luparm[LUSOL_IP_NONZEROS_L0],
                   LUSOL->luparm[LUSOL_IP_NONZEROS_U0],
                   (int) test, test/nnzero);

  test   = rnorm / (Amax*xnorm);
  fprintf(outunit, "\nAccuracy statistics:\n");
  fprintf(outunit, "%s with a factor tolerance of %g gave a relative error of %g\n",
                   LUSOL_pivotLabel(LUSOL), LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij], test);
  fprintf(outunit, "Amax = %g   bnorm = %g   rnorm = %g   xnorm = %g\n",
                   Amax, bnorm, rnorm, xnorm);
  fprintf(outunit, "\n");

  if (test <= 1.0e-8)
    fprintf(outunit, "The equations were solved with very high accuracy.\n");
  else if (test <= 1.0e-6)
    fprintf(outunit, "The equations were solved with reasonably good accuracy.\n");
  else {
    if (test <= 1.0e-4)
      fprintf(outunit, "Questionable accuracy; the LU factors may not be good enough.\n");
    else
      fprintf(outunit, "Poor accuracy; the LU factorization probably failed.\n");
    if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TRP)
      fprintf(outunit, "Try a smaller factor tolerance (current is %g).\n",
                       LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij]);
    else
      fprintf(outunit, "Try a smaller factor tolerance and/or TRP pivoting.\n");
  }

 /* Check if we should replace a column and resolve */
  if(replace > 0) {
    replace--;
    if(mod(replace, 2) == 1) {
      /* Randomly find a column and replace the data with the data in b */
      rnorm   = rand();
#ifdef YZHANG
      /* XXX: avoid integer overflow. */
      randcol = (int) (n * rnorm / (RAND_MAX+1.0)) + 1;
#else
      randcol = (int) (n * rnorm / (RAND_MAX+1)) + 1;
#endif      
#if 1
      MEMCLEAR(x, m+1);
      for(i = 1; i < m; i++)
        x[i] = Amax * rand() / RAND_MAX;
      inform = LUSOL_replaceColumn(LUSOL, randcol, x);
#else
      inform = LUSOL_replaceColumn(LUSOL, randcol, b);
#endif
    }
    else {
      /* Set the previously replaced column back and resolve */
      for (k = 1; k <= nnzero; k++) {
        i    = iA[k];
        j    = jA[k];
        if(j == randcol)
          x[i] = Aij[k];
      }
      inform = LUSOL_replaceColumn(LUSOL, randcol, x);
    }
    if(inform != LUSOL_INFORM_LUSUCCESS)
      fprintf(outunit, "Error:\n%s\n", LUSOL_informstr(LUSOL, inform));
    else
      goto Resolve;
  }


/* Free memory */
x900:
  if(!success)
    fprintf(outunit, "Insufficient memory or data file not found.\n");
  free(Aij);
  free(b);
  free(xexact);
  free(iA);
  free(jA);

  free(rhs);
  free(r);
  free(x);

/*  
//  LUSOL_dump(NULL, LUSOL);
// "D:\Program Files\CODE\Visual Studio Projects\LU\MatrixMarket\sherman5.mtx" "D:\Program Files\CODE\Visual Studio Projects\LU\MatrixMarket\sherman5_rhs1.mtx"
// A6805.txt b6805.txt
// A10009.txt b10009.txt
// fidap005.mtx fidap005_rhs1.mtx
// fidapm05.mtx fidapm05_rhs1.mtx
*/
  LUSOL_free(LUSOL);

/*     End of main program for Test of LUSOL */
  return(0);
}


