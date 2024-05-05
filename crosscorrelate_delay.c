#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "correlation.h"
#include "correl_nr.h"

/* 

   cross correlation offset

 */


int main(int argc,char **argv)
{
  COMP_PRECISION *t=NULL,*y1=NULL,*y2=NULL,*corr=NULL;
  COMP_PRECISION dt,corr_max,tcorr_max,rs,rsa;
  int i,j,n,nn,mlag;
  COMP_PRECISION rfac = 2;	/* subsample with 1/rfac  */

  if(argc < 3){
    fprintf(stderr,"%s file1 file2 [rfac, %g]\nreads in two x-y datafiles and computes cross-correlation as function of offset\n",
	    argv[0],rfac);
    exit(-1);
  }
  if(argc >= 4)
    sscanf(argv[3],"%lf",&rfac);
  if(rfac<1){
    fprintf(stderr,"%s: WARNING: sampling larger than unity?! %g\n",argv[0],rfac);
  }
 
  read_two_files_and_interpolate(&t,&y1, &y2,&n,rfac,(argv+1));
  dt=t[1]-t[0];
  compute_correl(&y1,&y2,&corr,n,&nn);
  corr_max=-1e20;
  
  mlag = nn/2;
  for(j=mlag,i= -mlag;j < nn;j++,i++){
    rs = corr[j]/(COMP_PRECISION)mlag;
    printf("%g %g\n",i*dt,rs);// negative lags
    rsa = fabs(rs);
    if(rsa > corr_max){
      corr_max = rsa;
      tcorr_max = i*dt;
    }
  }
  for(i=0;i <= mlag;i++){
    rs = corr[i]/(COMP_PRECISION)mlag;
    printf("%g %g\n",i*dt,rs);// positive lags
    rsa = fabs(rs);
    if(rsa > corr_max){
      corr_max = rsa;
      tcorr_max = i*dt;
    }
  }
  fprintf(stderr,"%s: max |r| correlation %g at t %g\n",argv[0],corr_max,tcorr_max);
  free(t);
  free(y1);
  free(y2);
  free(corr);
}

