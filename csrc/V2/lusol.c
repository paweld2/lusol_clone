
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   LUSOL routines from the Stanford Optimization Laboratory
   The parts included are:
    lusol1      Factor a given matrix A from scratch (lu1fac).
    lusol2      Heap-management routines for lu1fac.
    lusol6      Solve with the current LU factors.
    lusol7      Utilities for all update routines.
    lusol8      Replace a column (Bartels-Golub update).
   ------------------------------------------------------------------
   26 Apr 2002: TCP implemented using heap data structure.
   01 May 2002: lu1DCP implemented.
   07 May 2002: lu1mxc must put 0.0 at top of empty columns.
   09 May 2002: lu1mCP implements Markowitz with cols searched
                in heap order.
                Often faster (searching 20 or 40 cols) but more dense.
   11 Jun 2002: TRP implemented.
                lu1mRP implements Markowitz with Threshold Rook
                Pivoting.
                lu1mxc maintains max col elements  (was lu1max.)
                lu1mxr maintains max row elements.
   12 Jun 2002: lu1mCP seems too slow on big problems (e.g. memplus).
                Disabled it for the moment.  (Use lu1mar + TCP.)
   14 Dec 2002: TSP implemented.
                lu1mSP implements Markowitz with TSP.
   07 Mar 2003: character*1, character*2 changed to f90 form.
                Comments changed from * in column to ! in column 1.
                Comments kept within column 72 to avoid compiler
                warning.
   06 Mar 2004: Translation to C by Kjell Eikland with the addition
                of data wrappers, parametric constants, various 
                helper routines, and dynamic memory reallocation.
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
/* #include <varargs.h>  For UNIX 5 compatibility */
#include <string.h>
#include <float.h>
#include <math.h>
#include "lusol.h"
#include "myblas.h"


/* LUSOL Object creation and destruction */

void *clean_realloc(void *oldptr, int width, int newsize, int oldsize)
{
  newsize *= width;
  oldsize *= width;
  oldptr = realloc(oldptr, newsize);
  if(newsize > oldsize)
    memset((char *)oldptr+oldsize, '\0', newsize-oldsize);
  return(oldptr);
}

MYBOOL LUSOL_realloc_a(LUSOLrec *LUSOL, int newsize)
{
  int oldsize;

  if(newsize < 0)
    newsize = LUSOL->lena + max(abs(newsize), LUSOL_MINDELTA_a);

  oldsize = LUSOL->lena;
  LUSOL->lena = newsize;
  if(newsize > 0)
    newsize++;
  if(oldsize > 0)
    oldsize++;

  LUSOL->a    = (REAL *) clean_realloc(LUSOL->a,    sizeof(*(LUSOL->a)),
                                                    newsize, oldsize);
  LUSOL->indc = (int *)  clean_realloc(LUSOL->indc, sizeof(*(LUSOL->indc)),
                                                    newsize, oldsize);
  LUSOL->indr = (int *)  clean_realloc(LUSOL->indr, sizeof(*(LUSOL->indr)),
                                                    newsize, oldsize);
  if((newsize == 0) ||
     ((LUSOL->a != NULL) && (LUSOL->indc != NULL) && (LUSOL->indr != NULL)))
    return( TRUE );
  else
    return( FALSE );
}

MYBOOL LUSOL_expand_a(LUSOLrec *LUSOL, int *delta_lena, int *right_shift)
{
#ifdef StaticMemAlloc
  return( FALSE );
#else
  int LENA, NFREE, LFREE;
  
  /* Add expansion factor to avoid having to resize too often/too much;
     (exponential formula suggested by Michael A. Saunders) */
  LENA = LUSOL->lena;
  *delta_lena = (int) ((*delta_lena) * pow(1.5, fabs(*delta_lena)/LENA));

  /* Expand it! */
  if((*delta_lena <= 0) || !LUSOL_realloc_a(LUSOL, LENA+(*delta_lena)))
    return( FALSE );

  /* Make sure we return the actual memory increase of a */
  *delta_lena = LUSOL->lena-LENA;

  /* Shift the used memory area to the right */
  LFREE = *right_shift;
  NFREE = LFREE+*delta_lena;
  LENA  -= LFREE-1;
  MEMMOVE(LUSOL->a+NFREE,    LUSOL->a+LFREE,    LENA);
  MEMMOVE(LUSOL->indr+NFREE, LUSOL->indr+LFREE, LENA);
  MEMMOVE(LUSOL->indc+NFREE, LUSOL->indc+LFREE, LENA);

  /* Also return the new starting position for the used memory area of a */
  *right_shift  = NFREE;

  LUSOL->expanded_a++;
  return( TRUE );
#endif
}

MYBOOL LUSOL_realloc_r(LUSOLrec *LUSOL, int newsize)
{
  int oldsize;

  if(newsize < 0)
    newsize = LUSOL->maxm + max(abs(newsize), LUSOL_MINDELTA_rc);

  oldsize = LUSOL->maxm;
  LUSOL->maxm = newsize;
  if(newsize > 0)
    newsize++;
  if(oldsize > 0)
    oldsize++;

  LUSOL->lenr  = (int *) clean_realloc(LUSOL->lenr,  sizeof(*(LUSOL->lenr)), 
                                                     newsize, oldsize);
  LUSOL->ip    = (int *) clean_realloc(LUSOL->ip,    sizeof(*(LUSOL->ip)), 
                                                     newsize, oldsize);
  LUSOL->iqloc = (int *) clean_realloc(LUSOL->iqloc, sizeof(*(LUSOL->iqloc)), 
                                                     newsize, oldsize);
  LUSOL->ipinv = (int *) clean_realloc(LUSOL->ipinv, sizeof(*(LUSOL->ipinv)), 
                                                     newsize, oldsize);
  LUSOL->locr  = (int *) clean_realloc(LUSOL->locr,  sizeof(*(LUSOL->locr)), 
                                                     newsize, oldsize);

  if((newsize == 0) ||
     ((LUSOL->lenr != NULL) &&
      (LUSOL->ip != NULL) && (LUSOL->iqloc != NULL) &&
      (LUSOL->ipinv != NULL) && (LUSOL->locr != NULL))) {

#ifndef ClassicHamaxR
#ifdef AlwaysSeparateHamaxR
    if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TRP)
#endif
    {
      LUSOL->amaxr = (REAL *) clean_realloc(LUSOL->amaxr, sizeof(*(LUSOL->amaxr)),
                                                          newsize, oldsize);
      if((newsize > 0) && (LUSOL->amaxr == NULL))
        return( FALSE );
    }
#endif
    return( TRUE );
  }
  else
    return( FALSE );
}

MYBOOL LUSOL_realloc_c(LUSOLrec *LUSOL, int newsize)
{
  int oldsize;

  if(newsize < 0)
    newsize = LUSOL->maxn + max(abs(newsize), LUSOL_MINDELTA_rc);

  oldsize = LUSOL->maxn;
  LUSOL->maxn = newsize;
  if(newsize > 0)
    newsize++;
  if(oldsize > 0)
    oldsize++;

  LUSOL->lenc  = (int *)  clean_realloc(LUSOL->lenc,  sizeof(*(LUSOL->lenc)),
                                                      newsize, oldsize);
  LUSOL->iq    = (int *)  clean_realloc(LUSOL->iq,    sizeof(*(LUSOL->iq)),
                                                      newsize, oldsize);
  LUSOL->iploc = (int *)  clean_realloc(LUSOL->iploc, sizeof(*(LUSOL->iploc)),
                                                      newsize, oldsize);
  LUSOL->iqinv = (int *)  clean_realloc(LUSOL->iqinv, sizeof(*(LUSOL->iqinv)),
                                                      newsize, oldsize);
  LUSOL->locc  = (int *)  clean_realloc(LUSOL->locc,  sizeof(*(LUSOL->locc)),
                                                      newsize, oldsize);
  LUSOL->w     = (REAL *) clean_realloc(LUSOL->w,     sizeof(*(LUSOL->w)),
                                                      newsize, oldsize);
#ifdef LUSOLSafeFastUpdate
  LUSOL->vLU6L = (REAL *) clean_realloc(LUSOL->vLU6L, sizeof(*(LUSOL->vLU6L)),
                                                      newsize, oldsize);
#else
  LUSOL->vLU6L = LUSOL->w;
#endif

  if((newsize == 0) ||
     ((LUSOL->w != NULL) && (LUSOL->lenc != NULL) &&
      (LUSOL->iq != NULL) && (LUSOL->iploc != NULL) &&
      (LUSOL->iqinv != NULL) && (LUSOL->locc != NULL))) {

#ifndef ClassicHamaxR
    if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TCP) {
      LUSOL->Ha = (REAL *) clean_realloc(LUSOL->Ha,   sizeof(*(LUSOL->Ha)),
                                                      newsize, oldsize);
      LUSOL->Hj = (int *)  clean_realloc(LUSOL->Hj,   sizeof(*(LUSOL->Hj)),
                                                      newsize, oldsize);
      LUSOL->Hk = (int *)  clean_realloc(LUSOL->Hk,   sizeof(*(LUSOL->Hk)),
                                                      newsize, oldsize);
      if((newsize > 0) && 
         ((LUSOL->Ha == NULL) || (LUSOL->Hj == NULL) || (LUSOL->Hk == NULL)))
        return( FALSE );
    }
#endif
#ifndef ClassicdiagU
    if(LUSOL->luparm[LUSOL_IP_KEEPLU] == FALSE) {
      LUSOL->diagU = (REAL *) clean_realloc(LUSOL->diagU, sizeof(*(LUSOL->diagU)),
                                                          newsize, oldsize);
      if((newsize > 0) && (LUSOL->diagU == NULL))
        return( FALSE );
    }
#endif
      
    return( TRUE );
  }
  else
    return( FALSE );
}

LUSOLrec *LUSOL_create(FILE *outstream, int msgfil, int pivotmodel)
{
  LUSOLrec *newLU;

  newLU = (LUSOLrec *) calloc(1, sizeof(*newLU));
  if(newLU == NULL)
    return( newLU );

  newLU->luparm[LUSOL_IP_SCALAR_NZA]       = LUSOL_MULT_nz_a;
  newLU->outstream = outstream;
  newLU->luparm[LUSOL_IP_PRINTUNIT]        = msgfil;
  newLU->luparm[LUSOL_IP_PRINTLEVEL]       = LUSOL_MSG_SINGULARITY; 

  LUSOL_setpivotmodel(newLU, pivotmodel);

  newLU->parmlu[LUSOL_RP_GAMMA]            = LUSOL_DEFAULT_GAMMA;

  newLU->parmlu[LUSOL_RP_ZEROTOLERANCE]    = 3.0e-13;

  newLU->parmlu[LUSOL_RP_SMALLDIAG_U]      = /*3.7e-11;*/
  newLU->parmlu[LUSOL_RP_EPSDIAG_U]        = 3.7e-11;

  newLU->parmlu[LUSOL_RP_COMPSPACE_U]      = 3.0e+0;

  newLU->luparm[LUSOL_IP_MARKOWITZ_MAXCOL] = 5;
  newLU->parmlu[LUSOL_RP_MARKOWITZ_CONLY]  = 0.3e+0;
  newLU->parmlu[LUSOL_RP_MARKOWITZ_DENSE]  = 0.5e+0;

  newLU->luparm[LUSOL_IP_KEEPLU]           = TRUE;

  return( newLU );
}

MYBOOL LUSOL_sizeto(LUSOLrec *LUSOL, int init_r, int init_c, int init_a)
{
  if(LUSOL_realloc_a(LUSOL, init_a) &&
     LUSOL_realloc_r(LUSOL, init_r) &&
     LUSOL_realloc_c(LUSOL, init_c))
    return( TRUE );
  else
    return( FALSE );
}

char *LUSOL_pivotLabel(LUSOLrec *LUSOL)
{
  static /*const*/ char *pivotText[LUSOL_PIVOT_MAX+1] = 
  {"TPP", "TRP", "TCP", "TSP"};
  return(pivotText[LUSOL->luparm[LUSOL_IP_PIVOTTYPE]]);
}

void LUSOL_setpivotmodel(LUSOLrec *LUSOL, int pivotmodel)
{
  if((pivotmodel <= LUSOL_PIVOT_DEFAULT) || (pivotmodel > LUSOL_PIVOT_MAX))
    pivotmodel = LUSOL_PIVOT_TPP;
  LUSOL->luparm[LUSOL_IP_PIVOTTYPE]        = pivotmodel;

  /* Set default pivot tolerances
     (note that UPDATEMAX should always be <= FACTORMAX) */
  if(pivotmodel == LUSOL_PIVOT_TPP) {
#if 1
    LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij]    = 100.0;
    LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij]    =  10.0;
#else
    LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij]    =  10.0;
    LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij]    =  25.0;
#endif
  }
  if(pivotmodel == LUSOL_PIVOT_TRP) {
    LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij]    =  5.0;
    LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij]    =  5.0;
  }
  else {
    LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij]    = 10.0;
    LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij]    = 10.0;
  }
}

void LUSOL_tightenpivot(LUSOLrec *LUSOL)
{
  double newvalue;

  newvalue = sqrt(LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij]);
  LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij]    = newvalue;
  LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij]    = MIN(newvalue, 
                                                 LUSOL->parmlu[LUSOL_RP_UPDATEMAX_Lij]);
}


char *LUSOL_informstr(LUSOLrec *LUSOL, int inform)
{
  static /*const*/ char *informText[LUSOL_INFORM_MAX-LUSOL_INFORM_MIN+1] = 
  {"LUSOL_INFORM_RANKLOSS: Lost rank",
   "LUSOL_INFORM_LUSUCCESS: Success",
   "LUSOL_INFORM_LUSINGULAR: Singular A",
   "LUSOL_INFORM_LUUNSTABLE: Unstable factorization",
   "LUSOL_INFORM_ADIMERR: Row or column count exceeded",
   "LUSOL_INFORM_ADUPLICATE: Duplicate A matrix entry found",
   "LUSOL_INFORM_ANEEDMEM: Insufficient memory for factorization",
   "LUSOL_INFORM_FATALERR: Fatal internal error",
   "LUSOL_INFORM_NOPIVOT: Found no suitable pivot"};
  if(inform < LUSOL_INFORM_MIN || inform > LUSOL_INFORM_MIN)
    inform = LUSOL->luparm[LUSOL_IP_INFORM];
  return(informText[inform-LUSOL_INFORM_MIN]);
}

void LUSOL_clear(LUSOLrec *LUSOL, MYBOOL nzonly)
{
  int len;

  LUSOL->nelem = 0;
  if(!nzonly) {

   /* lena arrays */
    len = LUSOL->lena + LUSOL_ARRAYOFFSET;
    MEMCLEAR(LUSOL->a,    len);
    MEMCLEAR(LUSOL->indc, len);
    MEMCLEAR(LUSOL->indr, len);

   /* maxm arrays */
    len = LUSOL->maxm + LUSOL_ARRAYOFFSET;
    MEMCLEAR(LUSOL->lenr,  len);
    MEMCLEAR(LUSOL->ip,    len);
    MEMCLEAR(LUSOL->iqloc, len);
    MEMCLEAR(LUSOL->ipinv, len);
    MEMCLEAR(LUSOL->locr,  len);

#ifndef ClassicHamaxR
    if((LUSOL->amaxr != NULL)
#ifdef AlwaysSeparateHamaxR
       && (LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TRP)
#endif
      )
      MEMCLEAR(LUSOL->amaxr, len);
#endif

   /* maxn arrays */
    len = LUSOL->maxn + LUSOL_ARRAYOFFSET;
    MEMCLEAR(LUSOL->lenc,  len);
    MEMCLEAR(LUSOL->iq,    len);
    MEMCLEAR(LUSOL->iploc, len);
    MEMCLEAR(LUSOL->iqinv, len);
    MEMCLEAR(LUSOL->locc,  len);
    MEMCLEAR(LUSOL->w,     len);

    if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TCP) {
      MEMCLEAR(LUSOL->Ha,  len);
      MEMCLEAR(LUSOL->Hj,  len);
      MEMCLEAR(LUSOL->Hk,  len);
    }
#ifndef ClassicdiagU
    if(LUSOL->luparm[LUSOL_IP_KEEPLU] == FALSE) {
      MEMCLEAR(LUSOL->diagU, len);
    }
#endif
      
  }
}


MYBOOL LUSOL_assign(LUSOLrec *LUSOL, int iA[], int jA[], REAL Aij[], int nzcount, MYBOOL istriplet)
{
  int k, m, n, ij, kol;

  /* Adjust the size of the a structure */
  if(nzcount > (LUSOL->lena/LUSOL->luparm[LUSOL_IP_SCALAR_NZA]) &&
     !LUSOL_realloc_a(LUSOL, nzcount*LUSOL->luparm[LUSOL_IP_SCALAR_NZA]))
    return( FALSE );

  m = 0;
  n = 0;
  kol = 1;
  for(k = 1; k <= nzcount; k++) {
    /* First the row indicator */
    ij = iA[k];
    if(ij > m) {
      m = ij;
      if(m > LUSOL->maxm &&
         !LUSOL_realloc_r(LUSOL, -(m / LUSOL_MINDELTA_FACTOR + 1)))
        return( FALSE );
    }
    LUSOL->indc[k] = ij;

    /* Then the column indicator;
       Handle both triplet and column count formats */
    if(istriplet)
      ij = jA[k];
    else {
      if(k >= jA[kol])
        kol++;
      ij = kol;
    }
    if(ij > n) {
      n = ij;
      if(n > LUSOL->maxn &&
         !LUSOL_realloc_c(LUSOL, -(n / LUSOL_MINDELTA_FACTOR + 1)))
        return( FALSE );
    }
    LUSOL->indr[k] = ij;

    /* Lastly the matrix value itself */
    LUSOL->a[k] = Aij[k];
  }
  LUSOL->m = m;
  LUSOL->n = n;
  LUSOL->nelem = nzcount;
  return( TRUE );
}

int LUSOL_loadColumn(LUSOLrec *LUSOL, int iA[], int jA, REAL Aij[], int nzcount)
{
  int i, nz, k;

  nz = LUSOL->nelem;
  i = nz + nzcount;
  if(i > (LUSOL->lena/LUSOL->luparm[LUSOL_IP_SCALAR_NZA]) &&
     !LUSOL_realloc_a(LUSOL, i*LUSOL->luparm[LUSOL_IP_SCALAR_NZA]))
  return( -1 );

  k = 0;
  for(i = 1; i <= nzcount; i++) {
    if(Aij[i] == 0)
      continue;
    if(iA[i] <= 0 || iA[i] > LUSOL->m || 
       jA <= 0 || jA > LUSOL->n) {
      LUSOL_report(LUSOL, 0, "Variable index outside of set bounds (r:%d/%d, c:%d/%d)\n",
                             iA[i], LUSOL->m, jA, LUSOL->n);
      continue;
    }
    k++;
    nz++;
    LUSOL->a[nz]    = Aij[i];
    LUSOL->indc[nz] = iA[i];
    LUSOL->indr[nz] = jA;
  }
  LUSOL->nelem = nz;
  return( k );
}

void LUSOL_free(LUSOLrec *LUSOL)
{
  LUSOL_realloc_a(LUSOL, 0);
  LUSOL_realloc_r(LUSOL, 0);
  LUSOL_realloc_c(LUSOL, 0);
  free(LUSOL);
}

void LUSOL_report(LUSOLrec *LUSOL, int msglevel, char *format, ...)
{
  va_list ap;

  va_start(ap, format);
  if(LUSOL == NULL) {
    vfprintf(stderr, format, ap);
  }
  else if(msglevel >= 0  /*LUSOL->luparm[2]*/) {
    if(LUSOL->writelog != NULL) {
      char buff[255];

      vsprintf(buff, format, ap);
      LUSOL->writelog(LUSOL, LUSOL->loghandle, buff);
    }
    if(LUSOL->outstream != NULL) {
      vfprintf(LUSOL->outstream, format, ap);
      fflush(LUSOL->outstream);
    }
  }
  va_end(ap);
}

void LUSOL_timer(LUSOLrec *LUSOL, int timerid, char *text)
{
  LUSOL_report(LUSOL, -1, "TimerID %d at %s - %s\n", 
                          timerid, "", text);
}


int LUSOL_ftran(LUSOLrec *LUSOL, REAL b[], MYBOOL prepareupdate)
{
  int  inform;
  REAL *vector;

  if(prepareupdate)
    vector = LUSOL->vLU6L;
  else
    vector = LUSOL->w;
  MEMCPY(vector, b, LUSOL->n+1);
  LU6SOL(LUSOL, LUSOL_SOLVE_Aw_v, vector, b, &inform);
  return(inform);
}


int LUSOL_btran(LUSOLrec *LUSOL, REAL b[])
{
  int inform;

  MEMCPY(LUSOL->w, b, LUSOL->m+1);
  LU6SOL(LUSOL, LUSOL_SOLVE_Atv_w, b, LUSOL->w, &inform);
  return(inform);
}


int LUSOL_replaceColumn(LUSOLrec *LUSOL, int jcol, REAL v[])
{ 
  int  inform;
  REAL DIAG, VNORM;

  LU8RPC(LUSOL, LUSOL_UPDATE_OLDNONEMPTY, LUSOL_UPDATE_NEWNONEMPTY,
                jcol, v, NULL,
                &inform, &DIAG, &VNORM);
  LUSOL->replaced_c++;
  return( inform );
}

int LUSOL_findColumnPosition(LUSOLrec *LUSOL, int jcol)
/* The purpose oif this routine is to find the slack row/column that 
   was singular in the last unsuccessful column update; 
   zero is returned if the search was unsuccessful. 
   (Source is Michael A. Saunders in private communication to KE) */
{
  int j;

  for(j = LUSOL->m; j > 0; j--)
    if(LUSOL->iq[j] == jcol)
      break;
  if(j > 0)
    j = LUSOL->ip[j];
  return(j);
}

char relationChar(REAL left, REAL right)
{
  if(left > right)
    return('>');
  else if(left == right)
    return('=');
  else
    return('<');
}

/* Retrieve the core modules ordered by order of dependency */

#include "lusol2.c"      /* Heap management */
#include "lusol6a.c"     /* Singularity checking and equation solving */
#include "lusol1.c"      /* Factorization and core components */
#include "lusol7a.c"     /* Utility routines for updates */
#include "lusol8a.c"     /* Column update */


void LUSOL_dump(FILE *output, LUSOLrec *LUSOL)
{
  MYBOOL userfile = (MYBOOL) (output != NULL);

  if(!userfile)
    output = fopen("LUSOL.dbg", "w");

  blockWriteREAL(output, "a", LUSOL->a, 1, LUSOL->lena);
  blockWriteINT(output, "indc", LUSOL->indc, 1, LUSOL->lena);
  blockWriteINT(output, "indr", LUSOL->indr, 1, LUSOL->lena);

  blockWriteINT(output, "ip", LUSOL->ip, 1, LUSOL->m);
  blockWriteINT(output, "iq", LUSOL->iq, 1, LUSOL->n);
  blockWriteINT(output, "lenc", LUSOL->lenc, 1, LUSOL->n);
  blockWriteINT(output, "lenr", LUSOL->lenr, 1, LUSOL->m);

  blockWriteINT(output, "locc", LUSOL->locc, 1, LUSOL->n);
  blockWriteINT(output, "locr", LUSOL->locr, 1, LUSOL->m);
  blockWriteINT(output, "iploc", LUSOL->iploc, 1, LUSOL->n);
  blockWriteINT(output, "iqloc", LUSOL->iqloc, 1, LUSOL->m);

  blockWriteINT(output, "ipinv", LUSOL->ipinv, 1, LUSOL->m);
  blockWriteINT(output, "iqinv", LUSOL->iqinv, 1, LUSOL->n);

  if(!userfile)
    fclose(output);
}
