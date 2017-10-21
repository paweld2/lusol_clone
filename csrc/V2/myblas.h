
#include "commonlib.h"

#define BASE         1


#ifdef __cplusplus
extern "C" {
#endif

#define subvec(item) (item - 1)
/* int subvec( int item ); */

int submat( int nrowb, int row, int col );
int posmat( int nrowb, int row, int col );

void dscal (int n, REAL da, REAL *dx, int incx );
void dcopy ( int n, REAL *dx, int incx, REAL *dy, int incy );
void daxpy( int n, REAL da, REAL *dx, int incx, REAL *dy, int incy );
void dswap( int n, REAL *dx, int incx, REAL *dy, int incy );
REAL ddot ( int n, REAL *dx, int incx, REAL *dy, int incy );
void dload ( int n, REAL da, REAL *dx, int incx );
int idamax( int n, REAL *x, int is );
REAL dnormi( int n, REAL *x );

void randomseed(int *seeds);
void randomdens( int n, REAL *x, REAL r1, REAL r2, REAL densty, int *seeds);
void ddrand( int n, REAL *x, int incx, int *seeds );

#ifdef __cplusplus
}
#endif


