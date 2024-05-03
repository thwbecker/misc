#include "lmcurve.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* function to be fit */
double fit_erf_erf( double , const double * );
/* driver */
void fit_erf(int , double *, double *, double *, int ); /*  */
double erfinv( double );
/* model function: error function 

   d + b erf (x/c)
   
*/

double fit_erf_erf( double t, const double *p )
{
  return p[0] + p[1] * erf(t/p[2]);
}

void fit_erf(int m, double *t, double *y, double *par, int verbose)
{
  lm_control_struct control = lm_control_double;
  lm_status_struct status;
  const int n = 3;
  int i;
  control.verbosity = 0;
  if(verbose)
    printf( "Fitting ...\n" );

  /* now the call to lmfit */
  lmcurve( n, par, m, t, y, fit_erf_erf, &control, &status );
  if(par[2]<0){			/* make the third parameter always
				   positive */
    par[2] = -par[2];
    par[1] = -par[1];
  }
  if(verbose){
    printf( "Results:\n" );
    printf( "status after %d function evaluations:\n  %s\n",
	    status.nfev, lm_infmsg[status.outcome] );
  
    printf("obtained parameters:\n");
    for ( i = 0; i < n; ++i)
      printf("  par[%i] = %12g\n", i, par[i]);
    printf("obtained norm:\n  %12g\n", status.fnorm );
    
    printf("fitting data as follows:\n");
    for ( i = 0; i < m; ++i)
      printf( "  t[%2d]=%4g y=%6g fit=%10g residue=%12g\n",
      i, t[i], y[i], fit_erf_erf(t[i],par), y[i] - fit_erf_erf(t[i],par) );
    
  }
 }
    
/* 
   Function to calculate inverse error function.  Rational approximation 
   is used to generate an initial approximation, which is then improved to 
   full accuracy by two steps of Newton's method.  Code is a direct 
   translation of the erfinv m file in matlab version 2.0. 
   
   Author:  Gary L. Pavlis, Indiana University
   Date:  February 1996
*/
    
#define CENTRAL_RANGE 0.7
#define MAXDOUBLE DBL_MAX
    

double erfinv( double y)
{
  double x,z,num,dem; /*working variables */
  /* coefficients in rational expansion */
  double a[4]={ 0.886226899, -1.645349621,  0.914624893, -0.140543331};
  double b[4]={-2.118377725,  1.442710462, -0.329097515,  0.012229801};
  double c[4]={-1.970840454, -1.624906493,  3.429567803,  1.641345311};
  double d[2]={ 3.543889200,  1.637067800};
  if(fabs(y) > 1.0) 
    return (atof("NaN"));  /* This needs IEEE constant*/
  if(fabs(y) == 1.0) 
    return((copysign(1.0,y))*MAXDOUBLE); 

  if( fabs(y) <= CENTRAL_RANGE ) {
    z = y*y;
    num = (((a[3]*z + a[2])*z + a[1])*z + a[0]);
    dem = ((((b[3]*z + b[2])*z + b[1])*z +b[0])*z + 1.0);
    x = y*num/dem;
  }else if( (fabs(y) > CENTRAL_RANGE) && (fabs(y) < 1.0) ){
    z = sqrt(-log((1.0-fabs(y))/2.0));
    num = ((c[3]*z + c[2])*z + c[1])*z + c[0];
    dem = (d[1]*z + d[0])*z + 1.0;
    x = (copysign(1.0,y))*num/dem;
  }

  /* Two steps of Newton-Raphson correction */
   x = x - (erf(x) - y)/( (2.0/sqrt(M_PI))*exp(-x*x));
   x = x - (erf(x) - y)/( (2.0/sqrt(M_PI))*exp(-x*x));
  return(x);
}


