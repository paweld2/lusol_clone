
#ifndef REAL
#define REAL double
#endif
#ifndef MYBOOL
#define MYBOOL unsigned char
#endif
#ifndef my_mod
#define my_mod(n, m)           ((m) > 0 ? ((n) - (m) * (int) ((n)/(m))) : 0)
#endif


#ifdef F2C_INCLUDE
/* This is a compatibility data wrapper for regression testing */
typedef struct _LUSOLrec {

  /* General data - Don't need */

  /* Parameter storage arrays - Don't need */

  /* Arrays of length lena+1 */
  int    lena, nelem;
  int    *indc, *indr;
  REAL   *a;

  /* Arrays of length maxm+1 (row storage) */
  int    maxm, m;
  int    *lenr, *ip, *iqloc, *ipinv, *locr;

  /* Arrays of length maxn+1 (column storage) */
  int    maxn, n;
  int    *lenc, *iq, *iploc, *iqinv, *locc;
  REAL   *w;

} LUSOLrec;
#endif


void blockWriteINT(FILE *output, char *label, int *vector, int first, int last)
{
  int i, k = 0;

  fprintf(output, label);
  fprintf(output, "\n");
  for(i = first; i <= last; i++) {
    fprintf(output, " %5d", vector[i]);
    k++;
    if(my_mod(k, 12) == 0) {
      fprintf(output, "\n");
      k = 0;
    }
  }
  if(my_mod(k, 12) != 0)
    fprintf(output, "\n");
}

void blockWriteBOOL(FILE *output, char *label, MYBOOL *vector, int first, int last)
{
  int i, k = 0;

  fprintf(output, label);
  fprintf(output, "\n");
  for(i = first; i <= last; i++) {
    fprintf(output, " %1d", vector[i]);
    k++;
    if(my_mod(k, 36) == 0) {
      fprintf(output, "\n");
      k = 0;
    }
  }
  if(my_mod(k, 36) != 0)
    fprintf(output, "\n");
}
void blockWriteREAL(FILE *output, char *label, REAL *vector, int first, int last)
{
  int i, k = 0;

  fprintf(output, label);
  fprintf(output, "\n");
  for(i = first; i <= last; i++) {
    fprintf(output, " %18g", vector[i]);
    k++;
    if(my_mod(k, 4) == 0) {
      fprintf(output, "\n");
      k = 0;
    }
  }
  if(my_mod(k, 4) != 0)
    fprintf(output, "\n");
}

#ifdef F2C_INCLUDE
LUSOLrec *createLUSOL(int m, int n, int nelem, int lena, /* int *luparm, REAL *parmlu, */
                      REAL *a, int *indc, int *indr,
                      int *ip, int *iq, int *lenc, int *lenr,
                      int *locc, int *locr, int *iploc, int *iqloc,
                      int *ipinv, int *iqinv, REAL *w)
{
  LUSOLrec *LUSOL = calloc(1, sizeof(*LUSOL));

  LUSOL->m = m;
  LUSOL->n = n;
  LUSOL->nelem = nelem;
  LUSOL->lena = lena;

/*  LUSOL->luparm = &luparm-1; */
/*  LUSOL->parmlu = &parmlu-1; */

  LUSOL->a = a-1;
  LUSOL->indc = indc-1;
  LUSOL->indr = indr-1;

  LUSOL->ip = ip-1;
  LUSOL->iq = iq-1;
  LUSOL->lenc = lenc-1;
  LUSOL->lenr = lenr-1;

  LUSOL->locc = locc-1;
  LUSOL->locr = locr-1;
  LUSOL->iploc = iploc-1;
  LUSOL->iqloc = iqloc-1;

  LUSOL->ipinv = ipinv-1;
  LUSOL->iqinv = iqinv-1;

  LUSOL->w = w-1;

  return(LUSOL);

}
#endif

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
