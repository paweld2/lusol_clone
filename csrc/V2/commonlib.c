
#include <sys/types.h>

#ifdef INTEGERTIME
#include <time.h>
#else
#include <sys/timeb.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "commonlib.h"


/* Math functions */
int mod(int n, int d)
{
  return(n % d);
}

/* Array functions */
int findIndex(int target, int *attributes, int size, int offset)
{
  int focusPos, beginPos, endPos;
  int focusAttrib, beginAttrib, endAttrib;

 /* Set starting and ending index offsets */
  beginPos = offset;
  endPos = beginPos + size - 1;
  if(endPos < beginPos)
    return(-1);

 /* Do binary search logic based on a sorted attribute vector */
  focusPos = (beginPos + endPos) / 2;
  beginAttrib = attributes[beginPos];
  focusAttrib = attributes[focusPos];
  endAttrib   = attributes[endPos];

  while(endPos - beginPos > LINEARSEARCH) {
    if(beginAttrib == target) {
      focusAttrib = beginAttrib;
      endPos = beginPos;
    }
    else if(endAttrib == target) {
      focusAttrib = endAttrib;
      beginPos = endPos;
    }
    else if(focusAttrib < target) {
      beginPos = focusPos + 1;
      beginAttrib = attributes[beginPos];
      focusPos = (beginPos + endPos) / 2;
      focusAttrib = attributes[focusPos];
    }
    else if(focusAttrib > target) {
      endPos = focusPos - 1;
      endAttrib = attributes[endPos];
      focusPos = (beginPos + endPos) / 2;
      focusAttrib = attributes[focusPos];
    }
    else {
      beginPos = focusPos;
      endPos = focusPos;
    }
  }

 /* Do linear (unsorted) search logic */
  if(endPos - beginPos <= LINEARSEARCH) {

    /* CPU intensive loop; provide alternative evaluation models */
    if(!DOFASTMATH) {
      /* Do traditional indexed access */
      focusAttrib = attributes[beginPos];
      while((beginPos < endPos) && (focusAttrib < target)) {
        beginPos++;
        focusAttrib = attributes[beginPos];
      }
    }
    else {
      /* Do fast pointer arithmetic */
      int *attptr = attributes + beginPos;
      while((beginPos < endPos) && ((*attptr) < target)) {
        beginPos++;
        attptr++;
      }
      focusAttrib = (*attptr);
    }
  }

 /* Return the index if a match was found, or signal failure with a -1 */
  if(focusAttrib == target)             /* Found; return retrieval index      */
    return(beginPos);
  else if(focusAttrib > target)         /* Not found; last item               */
    return(-beginPos);
  else if(beginPos > offset+size-1)
    return(-(endPos+1));                /* Not found; end of list             */
  else
    return(-(beginPos+1));              /* Not found; intermediate point      */

}
int sortByREAL(int *item, REAL *weight, int size, int offset, MYBOOL unique)
{
  int i, ii, saveI;
  REAL saveW;

  for(i = 1; i < size; i++) {
    ii = i+offset-1;
    while ((ii >= 0) && (weight[ii] >= weight[ii+1])) {
      if(weight[ii] == weight[ii+1]) {
        if(unique)
          return(item[ii]);
      }
      else {
        saveI = item[ii];
        saveW = weight[ii];
        item[ii] = item[ii+1];
        weight[ii] = weight[ii+1];
        item[ii+1] = saveI;
        weight[ii+1] = saveW;
      }
      ii--;
    }
  }
  return(0);
}
int sortByINT(int *item, int *weight, int size, int offset, MYBOOL unique)
{
  int i, ii, saveI;
  int saveW;

  for(i = 1; i < size; i++) {
    ii = i+offset-1;
    while ((ii >= 0) && (weight[ii] >= weight[ii+1])) {
      if(weight[ii] == weight[ii+1]) {
        if(unique)
          return(item[ii]);
      }
      else {
        saveI = item[ii];
        saveW = weight[ii];
        item[ii] = item[ii+1];
        weight[ii] = weight[ii+1];
        item[ii+1] = saveI;
        weight[ii+1] = saveW;
      }
      ii--;
    }
  }
  return(0);
}


/* Time and message functions */
double timeNow()
{
#ifdef INTEGERTIME
  return((double)time(NULL));
#elif defined CLOCKTIME
  return((double)clock()/CLOCKS_PER_SEC /* CLK_TCK */);
#else
  struct timeb buf;
 
  ftime(&buf);
  return((double)buf.time+((double) buf.millitm)/1000.0);
#endif
}


/* Miscellaneous reporting functions */

/* List a vector of INT values for the given index range */
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

/* List a vector of MYBOOL values for the given index range */
void blockWriteBOOL(FILE *output, char *label, MYBOOL *vector, int first, int last, MYBOOL asRaw)
{
  int i, k = 0; 

  fprintf(output, label);
  fprintf(output, "\n");
  for(i = first; i <= last; i++) {
    if(asRaw)
      fprintf(output, " %1d", vector[i]);
    else
      fprintf(output, " %5s", my_boolstr(vector[i]));
    k++;
    if(my_mod(k, 36) == 0) {
      fprintf(output, "\n");
      k = 0;
    }
  }
  if(my_mod(k, 36) != 0)
    fprintf(output, "\n");
}

/* List a vector of REAL values for the given index range */
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


/* CONSOLE vector and matrix printing routines */
void printvec( int n, REAL *x, int modulo )
{
  int i;

  if (modulo <= 0) modulo = 5;
  for (i = 1; i<=n; i++) {
    if(mod(i, modulo) == 1) 
      printf("\n%2d:%12g", i, x[i]);
    else
      printf(" %2d:%12g", i, x[i]);
  }
  if(mod(i, modulo) != 0) printf("\n");
}


void printmatUT( int size, int n, REAL *U, int modulo)
{
   int i, ll;
   ll = 0;
   for(i = 1; i<=n; i++) {
     printvec(n-i+1, &U[ll], modulo);
     ll += size-i+1;
   }
}


void printmatSQ( int size, int n, REAL *X, int modulo)
{
   int i, ll;
   ll = 0;
   for(i = 1; i<=n; i++) {
     printvec(n, &X[ll], modulo);
     ll += size;
   }
}

