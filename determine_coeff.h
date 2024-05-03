#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "precision.h"
#include "splines.h"
#include "chebyshev.h"
#define MEMERROR {fprintf(stderr,"mem error\n");exit(-1);}


#define SPLINE_BASE 0
#define CHEBYSHEV_BASE 1
#define SPLINE_BASE_OUT 10
#define CHEBYSHEV_BASE_OUT 11


#define svdcmp svdcmp_
#define svbksb svbksb_
void svd_solver(COMP_PRECISION *, COMP_PRECISION *,
		COMP_PRECISION *,int , int );
void svd_driver(COMP_PRECISION *,int ,int ,
		COMP_PRECISION *,COMP_PRECISION *);

void svdcmp(COMP_PRECISION *,int *,int *,int *,
		   int *,COMP_PRECISION *,COMP_PRECISION *,
		   COMP_PRECISION *);
void svbksb(COMP_PRECISION *,COMP_PRECISION *,
		   COMP_PRECISION *,int *,int *,
		   int *,int *,COMP_PRECISION *,
		   COMP_PRECISION *,COMP_PRECISION *);
void determine_coeff(COMP_PRECISION *,int,COMP_PRECISION *,
		     COMP_PRECISION *, int ,
		     COMP_PRECISION , COMP_PRECISION ,
		     COMP_PRECISION , COMP_PRECISION ,
		     COMP_PRECISION *, 
		     COMP_PRECISION *,
		     COMP_PRECISION (*)(COMP_PRECISION,
					COMP_PRECISION,
					COMP_PRECISION *,int,
					COMP_PRECISION),
		     COMP_PRECISION (*)(int,int,int));
