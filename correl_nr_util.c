/*

  calculates the cross-correlation 
  of two timeseries
  using Numerical Recipes routines



*/
#include "correl_nr.h"
/* 

   compute correlation for different lags using FFT, and normalize by
   the two RMS values to actually get correlation

   corr nn on output
   

 */
void compute_correl(double *x, double *y, double **corr, int n, int *nn,int verbose)
{

  int lag;
  int i;
  double xm,ym,sum1,sum2,rmst,nn2,*xl,*yl;
  // max lag
  lag = n / 2;
  // zero padding
  *nn=pow(2,ceil(log2((double)(n+lag))));
  nn2 = *nn * 2;
  
  /* remove means */
  xm = ym = 0;
  sum1 = sum2 = 0.0;
  for(i=0;i<n;i++){
    xm += x[i];sum1 += x[i] * x[i];
    ym += y[i];sum2 += y[i] * y[i];
  }
  xm /= (double)n;
  ym /= (double)n;
  sum1 /= (double)n;
  sum2 /= (double)n;

  xl=(double *)calloc(*nn,sizeof(double));
  yl=(double *)calloc(*nn,sizeof(double));
  *corr=(double *)realloc(*corr,sizeof(double)*nn2);
  if(!xl || !yl || !(*corr))ME;
  
  for(i=0;i < n;i++){
    xl[i] = x[i] -= xm;
    yl[i] = y[i] -= ym;
  }
  if(verbose)
    fprintf(stderr,"compute_correl: using %i data pairs (removed means %g %g), next power of two padding for max lag of %i: %i\n",
	    n,xm,ym,lag,*nn);

  // calculate correlation, numrec style call
  correl((xl-1),(yl-1),(*nn),(*corr-1));
  /* normalize to get correlation */
  //fprintf(stderr,"normalizing by %g\n",rmst);
  rmst = sqrt(sum1*sum2)*(n);
  for(i=0;i<nn2;i++)
    *(*corr+i) /= rmst;
  free(xl);free(yl);
}
/*

  from here, numerical recipes routines

*/
void correl(double *data1,double *data2,int n,double *ans)
{
  unsigned long no2,i;
  double dum,*fft;
 
  
  fft=vector(1,n<<1);
  twofft(data1,data2,fft,ans,n);
  no2=n>>1;
  for (i=2;i<=n+2;i+=2) {
    ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;
    ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;
  }
  ans[2]=ans[n+1];
  realft(ans,n,-1);
  
  free_vector(fft,1,n<<1);
}

void twofft(double *data1,double *data2,double *fft1,double *fft2,
	    int n)
{
  unsigned long nn3,nn2,jj,j;
  double rep,rem,aip,aim;
  
  nn3=1+(nn2=2+n+n);
  for (j=1,jj=2;j<=n;j++,jj+=2) {
    fft1[jj-1]=data1[j];
    fft1[jj]=data2[j];
  }
  four1(fft1,n,1);
  fft2[1]=fft1[2];
  fft1[2]=fft2[2]=0.0;
  for (j=3;j<=n+1;j+=2) {
    rep=0.5*(fft1[j]+fft1[nn2-j]);
    rem=0.5*(fft1[j]-fft1[nn2-j]);
    aip=0.5*(fft1[j+1]+fft1[nn3-j]);
    aim=0.5*(fft1[j+1]-fft1[nn3-j]);
    fft1[j]=rep;
    fft1[j+1]=aim;
    fft1[nn2-j]=rep;
    fft1[nn3-j] = -aim;
    fft2[j]=aip;
    fft2[j+1] = -rem;
    fft2[nn2-j]=aip;
    fft2[nn3-j]=rem;
  }
}
void four1(double *data,int nn,int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
void realft(double *data,int n,int isign)
{
  unsigned long i,i1,i2,i3,i4,np3;
  double c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;
  
  theta=3.141592653589793/(double) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data,n>>1,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for (i=2;i<=(n>>2);i++) {
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    four1(data,n>>1,-1);
  }
}
/* allocate a double vector with subscript range v[nl..nh] */
double *vector(int nl,int nh)
{
  double *v;
  
  v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) {fprintf(stderr,"allocation failure in vector\n");exit(-1);}
  return v-nl+NR_END;
}

/* free a double vector allocated with vector() */
void free_vector(double *v,int nl,int nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}

