#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "correlation.h"



int read_two_files_and_interpolate(COMP_PRECISION **ti,COMP_PRECISION **y1, COMP_PRECISION **y2,
				   int *n,COMP_PRECISION rfac,char **filename)
{
  int i,j,use;
  COMP_PRECISION *t[2],*y[2],dtmin[2],tmin,tmax,dt,tloc,*yi[2];
  int n0[2];
  FILE *in;
  //fprintf(stderr,"%g %g\n",window_length,rfac);
  for(i=0;i < 2;i++){
    /* read in two time series */
    t[i]=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
    y[i]=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));

    n0[i]=0;
    in = fopen(filename[i],"r");
    if(!in){
      fprintf(stderr,"read_two_files_and_interpolate: can not open file %i %s\n",i+1,filename[i]);
      exit(-1);
    }
    dtmin[i] = 1e20;
    while(fscanf(in,"%lf %lf",(t[i]+n0[i]),(y[i]+n0[i]))==2){
      use = 1;
      if(finite(*(y[i]+n0[i]))){
	if(n0[i]){
	  dt = *(t[i]+n0[i]) - *(t[i]+n0[i]-1);
	  if(dt < 0){
	    fprintf(stderr,"read_two_files_and_interpolate: error file %s, entry %i, time (column one) has to be monotically increasing, ti %g %g ti-1 %g %g\n",
		    filename[i],n0[i]+1,*(t[i]+n0[i]),*(y[i]+n0[i]), *(t[i]+n0[i]-1), *(y[i]+n0[i]-1));
	    exit(-1);
	  }else if(dt==0){	/* if there's double entries, use old (not great) */
	    fprintf(stderr,"read_two_files_and_interpolate: WARNING file %s, entry %i, time is same, using old, ti %g %g ti-1 %g %g, and discarding this entry\n",
		    filename[i],n0[i]+1,*(t[i]+n0[i]),*(y[i]+n0[i]), *(t[i]+n0[i]-1), *(y[i]+n0[i]-1));
	    use = 0;
	  }else{
	    if(dt < dtmin[i])
	      dtmin[i]=dt;
	  }
	}
      }else{
	use = 0;
      }
      if(use){
	n0[i]++;
	t[i]=(COMP_PRECISION *)realloc(t[i],(n0[i]+1)*sizeof(COMP_PRECISION));
	y[i]=(COMP_PRECISION *)realloc(y[i],(n0[i]+1)*sizeof(COMP_PRECISION));
	if(!t[i] || !y[i])MEMERROR;
      }
    }
    fclose(in);
    fprintf(stderr,"read_two_files_and_interpolate: file %i: read %i values from %s, tmin/tmax: %g/%g dtmin: %g\n",
	    i+1,n0[i],filename[i],*(t[i]),*(t[i]+n0[i]-1),dtmin[i]);
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
  *n = ((int)((tmax-tmin)/dt))+1; /* timesteps after resampling */

  dt = (tmax-tmin)/((COMP_PRECISION)(*n)-1);
  
  fprintf(stderr,"read_two_files_and_interpolate: sampling n %i at dt %g tmin %g to tmax %g\n",
	  *n,dt,tmin,tmax);
 
  /* resample, i.e. interpolate */
  *y1=(COMP_PRECISION *)realloc(*y1,sizeof(COMP_PRECISION)*(*n));
  *y2=(COMP_PRECISION *)realloc(*y2,sizeof(COMP_PRECISION)*(*n));
  *ti=(COMP_PRECISION *)realloc(*ti,sizeof(COMP_PRECISION)*(*n));
  if(!(*y1) || !(*y2) || !(*ti))MEMERROR;
  yi[0] = *y1;yi[1] = *y2;

  tloc = tmin;
  for(i=0;i < 2;i++){
    for(j=0;j < *n;j++){
      if(i==0){
	*(*ti+j) = tloc;
	tloc += dt;
      }
      *(yi[i]+j) = interpolate(t[i],y[i],n0[i],*(*ti+j));
    }
  }
  tloc -= dt;
  if(fabs(tloc-tmax)>5e-6){	/* checck */
    fprintf(stderr,"read_two_files_and_interpolate: interpolation error, %g vs %g, %e\n",
	    tmax,tloc,fabs(tloc-tmax));
    exit(-1);
  }
  fprintf(stderr,"read_two_files_and_interpolate: sampling n %i at dt %g tmin %g to tmax %g\n",
	  *n,dt,tmin,tmax);
 
  for(i=0;i < 2;i++){
    free(t[i]);
    free(y[i]);
  }
  return 0;
}

COMP_PRECISION interpolate(COMP_PRECISION *x, COMP_PRECISION *y,int n, COMP_PRECISION x0)
{
  int j,k;
  COMP_PRECISION val,a,b;
  j=k=0;
  while(k < n && x[k] <= x0)
    k++;
  if(k==n){
    k--;j=k-1;
  }else if(k==0){
    j=0;k=1;
  }else{
    j=k-1;
  }
  // linearly interpolate
  a=(x[k] - x0)/(x[k] - x[j]);
  b=1.0-a;
  val  = a*y[j];// f[y[j]]
  val += b*y[k];// f[y[k]]=f[y[j+1]]
  return val;

}

COMP_PRECISION correlation(COMP_PRECISION *x, COMP_PRECISION *y, int n)
{
  int j;
  COMP_PRECISION yt,xt;
  COMP_PRECISION syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;
  
  for (j=0;j<n;j++) {
    ax += x[j];
    ay += y[j];
  }
  ax /= (COMP_PRECISION)n;
  ay /= (COMP_PRECISION)n;
  for (j=0;j<n;j++) {
    xt=x[j]-ax;
    yt=y[j]-ay;
    sxx += xt*xt;
    syy += yt*yt;
    sxy += xt*yt;
  }
  return sxy/sqrt(sxx*syy);
}
