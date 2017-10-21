
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>

#define BIGNUMBER    1.0e+30
#define TINYNUMBER   1.0e-4
#define MACHINEPREC  2.22e-16
#define MATHPREC     1.0e-16
#define ERRLIMIT     1.0e-6

#ifndef LINEARSEARCH
  #define LINEARSEARCH 5
#endif

#if 0
  #define INTEGERTIME
#endif

#ifndef MYBOOL
  #define MYBOOL       unsigned char
#endif

#ifndef REAL
  #define REAL         double
#endif

#ifndef my_mod
#if 0
  #define my_mod(n, m)   ((m) > 0 ? ((n) - (m) * (int) ((n)/(m))) : 0)
#endif
  #define my_mod(n, m)   ((n) % (m))
#endif

#ifndef my_boolstr
  #define my_boolstr(x)          (!(x) ? "FALSE" : "TRUE")
#endif

#ifndef NULL
  #define NULL 	       0 
#endif

#ifndef FALSE
  #define FALSE        0
  #define TRUE         1
#endif

#ifndef DOFASTMATH
  #define DOFASTMATH   TRUE
#endif


#ifndef CALLOC
#define CALLOC(ptr, nr)\
  if(!((void *) ptr = calloc((size_t)(nr), sizeof(*ptr))) && nr) {\
    printf("calloc of %d bytes failed on line %d of file %s\n",\
           (size_t) nr * sizeof(*ptr), __LINE__, __FILE__);\
  }
#endif

#ifndef MALLOC
#define MALLOC(ptr, nr)\
  if(!((void *) ptr = malloc((size_t)((size_t) (nr) * sizeof(*ptr)))) && nr) {\
    printf("malloc of %d bytes failed on line %d of file %s\n",\
           (size_t) nr * sizeof(*ptr), __LINE__, __FILE__);\
  }
#endif

#ifndef REALLOC
#define REALLOC(ptr, nr)\
  if(!((void *) ptr = realloc(ptr, (size_t)((size_t) (nr) * sizeof(*ptr)))) && nr) {\
    printf("realloc of %d bytes failed on line %d of file %s\n",\
           (size_t) nr * sizeof(*ptr), __LINE__, __FILE__);\
  }
#endif

#ifndef FREE
#define FREE(ptr)\
  if((void *) ptr != NULL) free(ptr);
#endif
  
#ifndef MEMCPY
#define MEMCPY(nptr, optr, nr)\
  memcpy((nptr), (optr), (size_t)((size_t)(nr) * sizeof(*(optr))))
#endif

#ifndef MEMMOVE
#define MEMMOVE(nptr, optr, nr)\
  memmove((nptr), (optr), (size_t)((size_t)(nr) * sizeof(*(optr))))
#endif

#ifndef MALLOCCCPY
#define MALLOCCPY(nptr, optr, nr)\
  {MALLOC(nptr, (size_t)(nr))\
   MEMCPY(nptr, optr, (size_t)(nr))}
#endif

#ifndef MEMCLEAR
#define MEMCLEAR(ptr, nr)\
  memset(ptr, '\0', (size_t)((size_t)(nr) * sizeof(*ptr)))
#endif

#define	ABS(x)	     ((x) < 0 ? -(x) : (x))
#define MIN(x, y)    ((x) < (y) ? (x) : (y))
#define MAX(x, y)    ((x) > (y) ? (x) : (y))
#define IF(t, x, y)  ((t) ? (x) : (y))
#define SIGN(x)      ((x) < 0 ? -1 : 1)


#ifdef __cplusplus
  extern "C" {
#endif


int mod( int n, int d );

int findIndex(int target, int *attributes, int size, int offset);
int sortByREAL(int *item, REAL *weight, int size, int offset, MYBOOL unique);
int sortByINT(int *item, int *weight, int size, int offset, MYBOOL unique);

double timeNow();

void blockWriteBOOL(FILE *output, char *label, MYBOOL *vector, int first, int last, MYBOOL asRaw);
void blockWriteINT(FILE *output, char *label, int *vector, int first, int last);
void blockWriteREAL(FILE *output, char *label, REAL *vector, int first, int last);

void printvec( int n, REAL *x, int modulo );
void printmatSQ( int size, int n, REAL *X, int modulo );
void printmatUT( int size, int n, REAL *U, int modulo );

#ifdef __cplusplus
  }
#endif
