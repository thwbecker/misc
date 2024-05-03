#include "lmcurve.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* function to be fit */
double fit_exp1( double , const double * );
double fit_exp2( double , const double * );
double fit_exp3( double , const double * );
/* driver */
double fit_exp_driver(int , double *, double *, double *, int, int,double (*)(double,const double *)); /*  */

