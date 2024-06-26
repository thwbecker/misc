#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/* 

   moving window cross correlation
   
*/

#define COMP_PRECISION double

#define MEMERROR {fprintf(stderr,"memory allocation error\n");exit(-1);}

COMP_PRECISION interpolate(COMP_PRECISION *, COMP_PRECISION *,int , COMP_PRECISION );
COMP_PRECISION correlation(COMP_PRECISION *,COMP_PRECISION *,int, int);
int read_two_files_and_interpolate(COMP_PRECISION **,COMP_PRECISION **, COMP_PRECISION **,
				   int *,COMP_PRECISION ,char **);
void find_max_from_nr_corr(COMP_PRECISION *,int ,COMP_PRECISION ,COMP_PRECISION *,COMP_PRECISION *,int);


void calc_lag_corr_pedestrian(COMP_PRECISION *,COMP_PRECISION *,
			      COMP_PRECISION **,int ,int *,int);
