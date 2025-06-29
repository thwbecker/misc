#include "momentmag.h"

/* get moment in Nm from magnitude */
double moment(double mag)
{
  const double three_over_two = 3./2.;
  return pow(10.,three_over_two * mag + 9.1);
}
/* get magnitude from moment in Nm */

double magnitude(double M0)
{
  const double two_over_three = 2./3.;
  return two_over_three *(log10(M0)-9.1);
}



void compute_histogram(float *mag,   /* magnitude vecors */
		       int n,	     /* number of magnitudes */
		       float *magmin, /* min and max magnitude range */
		       float *magmax, /* if 1e10 for min, or -1e10 for max, determine in routine */
		       int nbox, /* number of bins */
		       int cum, /* cumulative? */
		       float *be, /* right end of bin location */
		       int *nb,   /* number of counts in bin */
		       float *dx)  /* spacing between bins */

{
  int i,j,k,checksum,nlast;
  float xrange, *magcopy;
  if(nbox < 2){
    fprintf(stderr,"compute_histogram: error: nbox: %i\n",nbox);
    exit(-1);
  }
  /* fix ranges? */
  if(*magmin == 1e10)
    for(i=0;i<n;i++)
      if(mag[i]< *magmin)
	*magmin = mag[i];
  if(*magmax == -1e10)
    for(i=0;i < n;i++)
      if(mag[i]> *magmax)
	*magmax = mag[i];
  // boundaries
  xrange= *magmax - *magmin;
  *dx = xrange/(float)(nbox);
  //
  //
  // set up bins
  magcopy=(float *)malloc(n*sizeof(float));
  if(!magcopy){
    fprintf(stderr,"compute_histogram: mem error\n");
    exit(-1);
  }
  for(i=0;i < nbox;i++){
    be[i] = *magmin + (*dx) * (i+1);	// right boundary
    nb[i] =0;
  }

  /* sort */
  memcpy(magcopy,mag,sizeof(float)*n);
  qsort(magcopy,n,sizeof(float),(int(*)(const void *, const void *))cfunc);
  //
  // count
  //
  if(cum){    // cumulative statistics count
    nb[0] = n;
    i=0;j=1;
    while((i < n)&&(j < nbox)){
      while((magcopy[i] < be[j])&&(i<n))
	i++;
      nb[j] = n-i-2;
      if(nb[j] < 0)
	nb[j] = 0;
      j++;
    }
  }else{
    //
    // non-cumulative statistics 
    //
    j = k = checksum = 0;
    for(i=0;i < n;i++){
      if(magcopy[i] >= be[j]){
	nb[j] = k;
	checksum += k;
	j++;nlast=j;
	k = 0;
      }
      k++;
    }
    nb[nlast] = k;
    checksum += k;
    if(checksum != n){
      fprintf(stderr,"compute_histogram: count error: ninbox: %i n %i\n",j,n);
    }
  }
  free(magcopy);
}

void print_histogram(float *be,int *nb,int nbox,int n,float dx, int cum, FILE *out)
{
  int i;
  fprintf(out,"# N = %i mag: min %g max %g\n",n,be[0]-dx,be[nbox-1]);
  if(cum){
    fprintf(out,"# cumulative stats: M_0 n(M_0)/N n(M_0) M_l M_r\n");
  }else{
    fprintf(out,"# non-cumulative stats: M_0 n(M_0)/N n(M_0)\n");
  }
  for(i=0;i < nbox;i++)
    fprintf(out,"%.6e %.6e %12i %.6e %.6e\n",
	   moment(be[i]-dx/2),(float)nb[i]/(float)n,nb[i],
	   moment(be[i]-dx),moment(be[i]));
}

int cfunc(float *x1, float *x2)
{
  if(*x1 < *x2)
    return -1;
  else if(*x1 > *x2)
    return 1;
  else
    return 0;
}
