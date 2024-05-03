#include <stdio.h>
#include <math.h>
#include <string.h>


/* 

generate random numbers

given different compilation flags, will either produce 

UNIFORM defined:      uniformly distributed
GAUSSIAN defined:    Gaussian distributed
EXPONENTIAL defined: exponentially distributed
GAUSS_VEL defined:   gaussian distributed, but two values (east and north)
                     for sigma_east sigma_north corr_east_north input)


based on numerical recipes code

$Id: my_rand_routines.c,v 1.5 2005/02/23 03:04:21 becker Exp becker $


*/
/*  select the generator to use */
//#define GEN_TO_USE ran1
#define GEN_TO_USE ran2

#ifndef COMP_EPS
#define COMP_EPS 5.0e-15
#endif
#ifndef COMP_PRECISION
#define COMP_PRECISION double
#endif


COMP_PRECISION ran1(long *);
COMP_PRECISION ran2(long *);
COMP_PRECISION ran3(long *);
COMP_PRECISION gasdev(long *);
COMP_PRECISION expdev(long *);
COMP_PRECISION squared(COMP_PRECISION);
void ellipse_convert (COMP_PRECISION ,COMP_PRECISION ,COMP_PRECISION , COMP_PRECISION *,
		      COMP_PRECISION *,COMP_PRECISION *) ;
