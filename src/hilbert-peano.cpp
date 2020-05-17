#include <R.h>
#include "Rmath.h"

// https://www.rdocumentation.org/packages/HilbertCurve/versions/1.2.2/topics/HilbertCurve

extern "C" {

  void hilbertpeano(double *x, int *xlen, int *depth, int *setseed) {

    int n = xlen[0];
    int k = depth[0], i;

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();
    int nsur2 = n / 2;
    void d2xy ( int m, int d, int *x, int *y );
    double Rf_runif(double a, double b);
    double length = R_pow(4.0, (double)(k + 1)) - 1.0;
    int *ptsxbefore, *ptsybefore, *ptsxafter, *ptsyafter;
    ptsxbefore = new int[1];
    ptsybefore = new int[1];
    ptsxafter = new int[1];
    ptsyafter = new int[1];
    double D, dbefore, dafter, dd;
    for (i = 0; i < nsur2; i++) {
      x[i] = Rf_runif(0.0, length);
      dbefore = floor(x[i]);
      dafter = ceil(x[i]);
      dd = x[i] - dbefore;
      d2xy(k + 2, (int)dbefore, ptsxbefore, ptsybefore);
      d2xy(k + 2, (int)dafter, ptsxafter, ptsyafter);
      D = sqrt(R_pow(ptsxafter[0] - ptsxbefore[0], 2.0) + R_pow(ptsyafter[0] - ptsybefore[0], 2.0));
      x[i] = (double)ptsxbefore[0] + dd * (double)(ptsxafter[0] - ptsxbefore[0]) * D;
      x[i + nsur2] = (double)ptsybefore[0] + dd * (double)(ptsyafter[0] - ptsybefore[0]) * D;
    }
    
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] ptsxbefore;
    delete[] ptsybefore;
    delete[] ptsxafter;
    delete[] ptsyafter;
        

// We return
    return;

  }

/******************************************************************************/

void d2xy ( int m, int d, int *x, int *y )

/******************************************************************************/
/*
  Purpose:

    D2XY converts a 1D Hilbert coordinate to a 2D Cartesian coordinate.

  Modified:

    05 December 2015

  Parameters:

    Input, int M, the index of the Hilbert curve.
    The number of cells is N=2^M.
    0 < M.

    Input, int D, the Hilbert coordinate of the cell.
    0 <= D < N * N.

    Output, int *X, *Y, the Cartesian coordinates of the cell.
    0 <= *X, *Y < N.
*/
{
  void rot ( int n, int *x, int *y, int rx, int ry );
  int i4_power ( int i, int j );
  int n;
  int rx;
  int ry;
  int s;
  int t = d;

  n = i4_power ( 2, m );

  *x = 0;
  *y = 0;
  for ( s = 1; s < n; s = s * 2 )
  {
    rx = 1 & ( t / 2 );
    ry = 1 & ( t ^ rx );
    rot ( s, x, y, rx, ry );
    *x = *x + s * rx;
    *y = *y + s * ry;
    t = t / 4;
  }
  return;
}
/******************************************************************************/

int i4_power ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_POWER returns the value of I^J.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the base and the power.  J should be nonnegative.

    Output, int I4_POWER, the value of I^J.
*/
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      Rprintf("\n" );
      Rprintf("I4_POWER - Fatal error!\n" );
      Rprintf("  I^J requested, with I = 0 and J negative.\n" );
      return 1;
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      Rprintf("\n" );
      Rprintf("I4_POWER - Fatal error!\n" );
      Rprintf("  I^J requested, with I = 0 and J = 0.\n" );
      return 1;
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
/******************************************************************************/

void rot ( int n, int *x, int *y, int rx, int ry ) 

/******************************************************************************/
/*
  Purpose:

    ROT rotates and flips a quadrant appropriately.

  Modified:

    05 December 2015

  Parameters:

    Input, int N, the length of a side of the square.  N must be a power of 2.

    Input/output, int *X, *Y, the old and the new coordinates.

    Input, int RX, RY, ???
*/
{
  int t;

  if ( ry == 0 )
  {
/*
  Reflect.
*/
    if ( rx == 1 )
    {
      *x = n - 1 - *x;
      *y = n - 1 - *y;
    }
/*
  Flip.
*/
     t = *x;
    *x = *y;
    *y =  t;
  }
  return;
}
/******************************************************************************/


int xy2d ( int m, int x, int y )

/******************************************************************************/
/*
  Purpose:

    XY2D converts a 2D Cartesian coordinate to a 1D Hilbert coordinate.

  Discussion:

    It is assumed that a square has been divided into an NxN array of cells,
    where N is a power of 2.

    Cell (0,0) is in the lower left corner, and (N-1,N-1) in the upper 
    right corner.

  Modified:

    05 December 2015

  Parameters:

    Input, int M, the index of the Hilbert curve.
    The number of cells is N=2^M.
    0 < M.

    Input, int X, Y, the Cartesian coordinates of a cell.
    0 <= X, Y < N.

    Output, int XY2D, the Hilbert coordinate of the cell.
    0 <= D < N * N.
*/
{
  int d = 0;
  int n;
  int rx;
  int ry;
  int s;

  n = i4_power ( 2, m );

  for ( s = n / 2; s > 0; s = s / 2 )
  {
    rx = ( x & s ) > 0;
    ry = ( y & s ) > 0;
    d = d + s * s * ( ( 3 * rx ) ^ ry );
    rot ( s, &x, &y, rx, ry );
  }
  return d;
}

  
}
