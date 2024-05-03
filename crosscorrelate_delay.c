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
  COMP_PRECISION tmin,tmax,*t[2],*y[2],*yi[2],dt,dtmin[2],tloc,rfac=1,*corr,corr_max,tcorr_max;
  FILE *in;
  int i,n0[2],n,j,mlag,nn,use;

  
  if(argc < 3){
    fprintf(stderr,"%s file1 file2 [rfac, %g]\nreads in two x-y datafiles and computes cross-correlation as function of offset\n",
	    argv[0],rfac);
    exit(-1);
  }
  if(argc>=4)
    sscanf(argv[3],"%lf",&rfac);
  for(i=0;i<2;i++){
    t[i]=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
    y[i]=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));

    n0[i]=0;
    in=fopen(argv[1+i],"r");
    if(!in){
      fprintf(stderr,"%s: can not open file %i %s\n",argv[0],i+1,argv[1+i]);
      exit(-1);
    }
    dtmin[i] = 1e20;
    while(fscanf(in,"%lf %lf",(t[i]+n0[i]),(y[i]+n0[i]))==2){
      use = 1;
      if(n0[i]){
	dt = *(t[i]+n0[i]) - *(t[i]+n0[i]-1);
	if(dt<0){
	  fprintf(stderr,"%s: error file %s, entry %i, time has to be monotically increasing, ti %g %g ti-1 %g %g\n",
		  argv[0],argv[i+1],n0[i]+1,*(t[i]+n0[i]),*(y[i]+n0[i]), *(t[i]+n0[i]-1), *(y[i]+n0[i]-1));
	  exit(-1);
	}else if(dt==0){
	  fprintf(stderr,"%s: WARNING file %s, entry %i, time is same, using old, ti %g %g ti-1 %g %g\n",
		  argv[0],argv[i+1],n0[i]+1,*(t[i]+n0[i]),*(y[i]+n0[i]), *(t[i]+n0[i]-1), *(y[i]+n0[i]-1));
	  use = 0;
	}else{
	  if(dt<dtmin[i])
	    dtmin[i]=dt;
	}
      }
      if(use){
	n0[i]++;
	t[i]=(COMP_PRECISION *)realloc(t[i],(n0[i]+1)*sizeof(COMP_PRECISION));
	y[i]=(COMP_PRECISION *)realloc(y[i],(n0[i]+1)*sizeof(COMP_PRECISION));
	if(!t[i] || !y[i])MEMERROR;
      }
    }
    fclose(in);
    fprintf(stderr,"%s: read %i values from %s, tmin/tmax: %g/%g dtmin: %g\n",argv[0],n0[i],argv[i+1],*(t[i]),*(t[i]+n0[i]-1),dtmin[i]);
  }
  if(*(t[0]) < *(t[1]))
    tmin= *(t[1]);
  else
    tmin= *(t[0]);
  if(*(t[0]+n0[0]-1) < *(t[1]+n0[1]-1))
    tmax = *(t[0]+n0[0]-1);
  else
    tmax = *(t[1]+n0[1]-1);
  if(dtmin[0] < dtmin[1])
    dt = dtmin[0]/rfac;
  else
    dt = dtmin[1]/rfac;
  n = ((int)((tmax-tmin)/dt))+1;

  dt = (tmax-tmin)/(n-1);
  
  fprintf(stderr,"%s: sampling n %i at dt %g  tmin %g to tmax %g\n",
	  argv[0],n,dt,tmin,tmax);
  /* resample */
  yi[0]=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n);
  yi[1]=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n);
  corr=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
  if(!yi[0] || !yi[1])MEMERROR;
  for(i=0;i < 2;i++){
    tloc = tmin;
    for(j=0;j < n;j++){
      *(yi[i]+j) = interpolate(t[i],y[i],n0[i],tloc);
      tloc += dt;
    }
  }
  tloc -= dt;
  if(fabs(tloc-tmax)>5e-6){	/* checck */
    fprintf(stderr,"%s: interpolation error, %g vs %g, %e\n",argv[0],tmax,tloc,fabs(tloc-tmax));
    exit(-1);
  }
  compute_correl(&(yi[0]),&(yi[1]),&corr,n,&nn);
  corr_max=-1e20;
  
  mlag = nn/2;
  for(j=mlag,i=-mlag;j<nn;j++,i++){
    printf("%g %g\n",i*dt,corr[j]);// negative lags
    if(corr[j]>corr_max){
      corr_max = corr[j];
      tcorr_max = i*dt;
    }
  }
  for(i=0;i<=mlag;i++){
    printf("%g %g\n",i*dt,corr[i]);// positive lags
    if(corr[i] > corr_max){
      corr_max = corr[i];
      tcorr_max = i*dt;
    }
  }
  fprintf(stderr,"%s: max correlation %g at t %g\n",argv[0],corr_max,tcorr_max);
  for(i=0;i<2;i++){
    free(t[i]);
    free(y[i]);
    free(yi[i]);
  }
  free(corr);
}

