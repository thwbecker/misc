
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define CPREC float
#define BOOLEAN unsigned short
#define FCPREC "%f"
#define FCPREC_2 "%f %f"

#define TRUE 1
#define FALSE 0


#define LAPACK_LLS sgels_  /* LAPACK least squares solver */

#ifndef MEMERROR
#define MEMERROR(x) {fprintf(stderr,"memory allocation error: %s\n",x);exit(-1);}
#endif
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))


extern void sgels_(char *,int *, int *, int *, float *,int *,
		   float *,int *,float *,int *,int *);

void solver_ab_lls(CPREC *, int , int , CPREC *);
void fit2dcircle(CPREC *, BOOLEAN *, int , int , CPREC *, CPREC *);
int d_compare(const void *,const void *);
