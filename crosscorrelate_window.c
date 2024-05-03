#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "correlation.h"
/* 

   moving window cross correlation

 */


int main(int argc,char **argv)
{
  COMP_PRECISION *t=NULL,*y1=NULL,*y2=NULL;
  int i,n,nr,nr2,il,imin,imax;
  COMP_PRECISION window_length = 0.2;	/* window length  */
  COMP_PRECISION rfac = 2;	/* subsample with 1/rfac  */
  
  if(argc < 3){
    fprintf(stderr,"%s file1 file2 [wl, %g] [rfac, %g]\nreads in two x-y datafiles and computes moving window (fraction wl) cross correlation\n\tresampling at dtmin/rfac\n",
	    argv[0],window_length,rfac);
    fprintf(stderr,"\tif wl < 0, will use actual time interval\n");
    
    exit(-1);
  }
  if(argc>=4)			/* windows length for correlation */
    sscanf(argv[3],"%lf",&window_length);
  if(argc>=5)			/* resampling factor */
    sscanf(argv[4],"%lf",&rfac);
  if(rfac>1){
    fprintf(stderr,"%s: WARNING: sampling larger than unity?! %g\n",argv[0],rfac);
  }
  read_two_files_and_interpolate(&t,&y1, &y2,&n,rfac,(argv+1));
  if(window_length < 0){		/* provided as actual time */
    fprintf(stderr,"%s: window length desired in actual time: %g\n",argv[0],-window_length);
    window_length = -window_length/(t[n-1]-t[0]);
  }
  if((window_length < 0)||(window_length>1)){
    fprintf(stderr,"%s: window length error %g (0....1)\n",argv[0],window_length);
    exit(-1);
  }

  nr = (int) ((COMP_PRECISION)n*window_length+0.5);
  if(nr%2!=0)
    nr++;

 
  /* output */
  nr2 = nr/2;
  il = 0;
  imin = il+nr2;imax = n - nr2;
  for(i=imin;i < imax;i++,il++){
    /* time correlation y1_int y2_int */
    printf("%20.7e %20.15f\t %g %g\n",
	   t[i],correlation((y1+i-nr2),(y2+i-nr2),nr),y1[i],y2[i]);
  }
  free(t);
  free(y1);
  free(y2);
}

