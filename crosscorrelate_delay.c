#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "correlation.h"
#include "correl_nr.h"

/* 

   cross correlation as a function of offset

 */


int main(int argc,char **argv)
{
  COMP_PRECISION *t=NULL,*y1=NULL,*y2=NULL,*corr=NULL;
  COMP_PRECISION dt,corr_max,tcorr_max;
  int i,j,n,nn,mlag;
  COMP_PRECISION rfac = 1;	/* subsample with 1/rfac  */
  int mode = 2;			/* 1: NR 2: mine 3: mine, pearson */
  if(argc < 3){
    fprintf(stderr,"%s file1 file2 [rfac, %g] [mode, %i]\nreads in two x-y datafiles and computes cross-correlation as function of offset\n",
	    argv[0],rfac,mode);
    exit(-1);
  }
  if(argc >= 4)
    sscanf(argv[3],"%lf",&rfac);
  if(rfac<1){
    fprintf(stderr,"%s: WARNING: sampling larger than unity?! %g\n",argv[0],rfac);
  }
  if(argc >= 5)
    sscanf(argv[4],"%i",&mode);
 
  read_two_files_and_interpolate(&t,&y1, &y2,&n,rfac,(argv+1));
  dt=t[1]-t[0];
  if(mode == 1)
    compute_correl(y1,y2,&corr,n,&nn,1); /* sorted differently on output */
  else
    calc_lag_corr_pedestrian(y1,y2,&corr,n,&nn,mode); /* simple sorting on output */
 
  find_max_from_nr_corr(corr,nn,dt,&tcorr_max,&corr_max,mode);
  fprintf(stderr,"%s: mode %i max |r| correlation %g at t %g\n",argv[0],mode,corr_max,tcorr_max);
  mlag = nn/2;
  if(mode==1){
    for(j=mlag,i= -mlag;j < nn;j++,i++){
      printf("%g %g\n",i*dt,corr[j]);// negative lags
    }
    for(i=0;i <= mlag;i++){
      printf("%g %g\n",i*dt,corr[i]);// positive lags
    }
  }else{
    for(j=-mlag,i=0;i < nn;i++,j++)
      printf("%g %g\n",j*dt,corr[i]);
  }
  free(t);
  free(y1);
  free(y2);
  free(corr);
}

