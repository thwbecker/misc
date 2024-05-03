#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* 

   moving window cross correlation
   
*/

#define COMP_PRECISION double
#define MEMERROR {fprintf(stderr,"memory allocation error\n");exit(-1);}

COMP_PRECISION interpolate(COMP_PRECISION *, COMP_PRECISION *,int , COMP_PRECISION );
COMP_PRECISION correlation(COMP_PRECISION *,COMP_PRECISION *,int );
int read_two_files_and_interpolate(COMP_PRECISION **,COMP_PRECISION **, COMP_PRECISION **,
				   int *,COMP_PRECISION ,char **);


