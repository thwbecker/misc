#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*

  implementation of the numerical recipes periodogram
  estimator based on uneven sampling using the Lomb 
  method

  input is time_i x_i from stdin, unformatted

  output is frequency power/max_power power error

  call with possible three arguments: 

  oversampling fraction_of_nyquist method

  method: 1 slow 2 fast periodogram 
          3 multitaper from F90 MTSPEC package
	  4 mutlitaper from JTAP Lees and Park MT


	  

  $Id: period.c,v 1.1 2002-12-15 14:47:21-08 tbecker Exp tbecker $

*/
/*

  macros

 */
#define COMP_PRECISION double	/* should be double for mt library */
#define TWOPID 6.2831853071795865
static double sqrargsdf;
#define SQUARE(x) (((sqrargsdf=(x))== 0.0) ? (0.0) : (sqrargsdf*sqrargsdf))
#define ABS(x) (((x)>=0)?((x)):(-(x)))
#define MEMERROR {fprintf(stderr,"%s: memory error\n",argv[0]);exit(-1);}
#define NR_END 1
#define FREE_ARG char*
#define MOD(a,b) 	while(a >= b) a -= b;
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

COMP_PRECISION lin_interpolate (COMP_PRECISION *, COMP_PRECISION *, int  , COMP_PRECISION);

#define MT_SPEC_FUNC mtspec_wrapper_
typedef struct {double r,i;} Fcomplex; /* should be float according to
					  web, but double makes more
					  sense */

void MT_SPEC_FUNC(int *,double *,double *,double *,int *,int *,double *,double *,unsigned char *,
		  int *,int *,Fcomplex *, double *, double *, double *, double *,
		  int *, double *, float *, int *);

void jlpmtm_do_mtap_spec(double *, int, int, int, double, int, double,
			double *, double *, double *, int); 
int jlpmtm_get_pow_2(int);

// for fasper: number of interpolation points per 1/4 cycle
#define MACC 4
/*

  function declarations

*/
void period(COMP_PRECISION *,COMP_PRECISION *, int ,
	    COMP_PRECISION , COMP_PRECISION ,COMP_PRECISION *,
	    COMP_PRECISION *,int ,int *,int *,COMP_PRECISION *);
void avevar(COMP_PRECISION *, long ,
	    COMP_PRECISION *,COMP_PRECISION *);
void nrerror(char *);
COMP_PRECISION *vector(long ,long );
void free_vector(COMP_PRECISION *,long ,long );
void spread(COMP_PRECISION ,COMP_PRECISION *,int ,
	    COMP_PRECISION ,int );
void four1(COMP_PRECISION *,int ,int );
void realft(COMP_PRECISION *,int ,int );
void fasper(COMP_PRECISION *,COMP_PRECISION *,int ,
       COMP_PRECISION ,COMP_PRECISION ,COMP_PRECISION *,
       COMP_PRECISION *,int ,int *,
       int *,COMP_PRECISION *);
/*

  start main

*/
int main(int argc, char **argv)
{
  int n,i,nout,jmax,nsample,nf;
  COMP_PRECISION *x,*y,fx,fc,fh,np,xmin,xmax,dx,dxmin,
    *px,prob,mean,max_spec;
  COMP_PRECISION ofac = 8.;	  /* for periodogram: oversampling for MT: resampling*/
  COMP_PRECISION hifac = 0.25;	  /* for periodogram hifac, MT: time-bandwith product, use 3.5 or 3 */

  /* default */
  int method = 2;		/* default is periodogram */
  /* for MTSPEC */
  double *wt,*err,*se,*sk,*fstat,*py,*freq,*spec;
  Fcomplex *yk;
  /* for LMPTM */
  double *dof,nyquist,df;
  int increase = 1;
  int klen;
  /* for MTSPEC */

  unsigned char verb[1] = "n";
  int qispec = 0;
  int adapt = 1;
  int nodemean = 1;		/* we already demean */
  int rshape = 0;
  float fcrit = 0.95;		/* weirdly, a real*4 */
  /* for both MT methods */
  int kspec = 5;		  /* number of tapers */
  /*

    parameters are 

    ofac:  oversampling: typically 4 or larger
    hifac: hifac*average Nyquist frequency is the highest 
           output frequency
    method: 1: old, slow 2: new, fast 3: MT MTSPEC 4: LMPTM
    
   */
  // oversampling
  if(argc >= 2)
    sscanf(argv[1],"%lf",&ofac);
  if(argc >= 3)
    sscanf(argv[2],"%lf",&hifac);
  if(argc >= 4)
    sscanf(argv[3],"%i",&method);
  /*

    reading data

  */ 
  n=0;
  x=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
  y=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
  fprintf(stderr,"%s: reading x y pairs from stdin\n",
	  argv[0]);
  mean=0.0;
  while(fscanf(stdin,"%lf %lf",(x+n),(y+n))==2){
    mean += y[n];
    n++;
    if((x=(COMP_PRECISION *)
	realloc(x,sizeof(COMP_PRECISION)*(n+1)))==NULL)
      MEMERROR;
    if((y=(COMP_PRECISION *)
	realloc(y,sizeof(COMP_PRECISION)*(n+1)))==NULL)
      MEMERROR;
  }
  xmin=1e20;xmax=-1e20;dxmin=1e20;
  for(i=0;i<n;i++){
    if(i>0){
      dx = x[i]-x[i-1];
      if(dx <= 0){
	fprintf(stderr,"%s: error, x values need to be sorted ascendingly\n",argv[0]);
	fprintf(stderr,"%s: line %i %g line %i %g dx %g\n",argv[0],i,x[i-1],i+1,x[i],dx);
	exit(-1);
      }
      if(dx < dxmin)
	dxmin = dx;
    }
    if(x[i]<xmin)
      xmin = x[i];
    if(x[i]>xmax)
      xmax = x[i];
  }
  if(n > 3){
    mean /= (COMP_PRECISION)n;
    fprintf(stderr,"%s: read %i data points, x: %g - %g - %g, removing y mean of %g\n",
	    argv[0],n,xmin,dxmin,xmax,mean);
    for(i=0;i<n;i++)
      y[i] -= mean;
    
    if(method <= 2){
      
      //
      // reset start of x to zero time
      fx=x[0];
      for(i=0;i<n;i++)
	x[i] -= fx;
      fprintf(stderr,"%s: periodogram oversampling factor: %g method: %i\n",argv[0],ofac,method);
      //
      // average nyquist frequency
      fc=(COMP_PRECISION)n/(2.0*x[n-1]);
      fh = hifac * fc;
      fprintf(stderr,"%s: avg nyquist f_c: %g f_high/f_c: %g\n",
	      argv[0],fc,fh/fc);
      // size of work array
      if(method == 1){
	np=10+(int)((ofac*hifac)/2.0*(COMP_PRECISION)n);
      }else if(method == 2){
	np=10+(int)((ofac*hifac)/2.0*(COMP_PRECISION)n)*ofac*8;
      }else{
	fprintf(stderr,"%s: method %i undefined\n",argv[0],method);
	exit(-1);
      }
      if((px=(COMP_PRECISION *)calloc(np,sizeof(COMP_PRECISION)))==NULL)
	MEMERROR;
      if((py=(COMP_PRECISION *)calloc(np,sizeof(COMP_PRECISION)))==NULL)
	MEMERROR;
      /*
	
	calculate Lomb periodogram
	
	on output, px will have the frequencies and py the value
	of the Lomb periodogram
	
	nout: the number of assigned frequencies
	py[jmax] is the maximum element of py[] and 
	prob, the significance of that max against random nois
	
	
      */
      if(method == 1){
	// first method
	period(x-1,y-1,n,ofac,hifac,px-1,py-1,np,&nout,&jmax,&prob);
      }else if(method == 2){
	// second, fast method
	fasper(x-1,y-1,n,ofac,hifac,px-1,py-1,np,&nout,&jmax,&prob);
      }
      
      jmax--;
      //fprintf(stdout,"# max: frequency: %g power: %g probability: %g%%\n", py[jmax],px[jmax],(1.0-prob)*100.);
      for(i=0;i<nout;i++)
	// output is frequency power/max_power power  NaN NaN
	fprintf(stdout,"%20.15e %20.15e %20.15e NaN NaN\n",
		px[i],py[i]/py[jmax],py[i]);
      free(px);free(py);
    }else{
    
      /* resample */
      dx = dxmin/ofac;
      nsample = (int)((xmax-xmin)/dx+0.5) + 1;
      if(nsample%2!=0)
	nsample++;		/* make sure it's even */
      dx = (xmax - xmin)/(nsample-1);
      fprintf(stderr,"%s: MT resampling factor: %g dx %g, tbw: %g\n",argv[0],ofac,dx,hifac);
      if((hifac<3)||(hifac> 3.5))
	fprintf(stderr,"%s: WARNING: strange value for tbw\n",argv[0]);
      
      if((py=(COMP_PRECISION *)calloc(nsample,sizeof(COMP_PRECISION)))==NULL)
	MEMERROR;
      for(i=0;i < nsample;i++){
	py[i] = lin_interpolate (x,y,n,xmin + i*(COMP_PRECISION)dx);
	//fprintf(stdout,"%g %g\n",xmin + i*(COMP_PRECISION)dx,py[i]);
      }

      if(method == 3){
	nf = nsample/2+1;
	if((freq=(double *)calloc(nf,sizeof(double)))==NULL)
	  MEMERROR;
	if((spec=(double *)calloc(nf,sizeof(double)))==NULL)
	  MEMERROR;
	
	/* optional */
	if((yk=(Fcomplex *)calloc(nsample * kspec,sizeof(Fcomplex)))==NULL)
	  MEMERROR;
	if((wt=(double *)calloc(nf*kspec,sizeof(double)))==NULL)
	  MEMERROR;
	if((err=(double *)calloc(nf*2,sizeof(double)))==NULL)
	  MEMERROR;
	if((se=(double *)calloc(nf,sizeof(double)))==NULL)
	  MEMERROR;
	if((sk=(double *)calloc(nf*kspec,sizeof(double)))==NULL)
	  MEMERROR;
	if((fstat=(double *)calloc(nf,sizeof(double)))==NULL)
	  MEMERROR;
	
	fprintf(stderr,"%s: MT SPEC kspec: %i\n",argv[0],kspec);
	/* multitaper */
	MT_SPEC_FUNC(&nsample,&dx,py,&hifac,&kspec,&nf,freq,spec,verb,
		     &qispec,&adapt,yk, wt, err, se, sk,
		     &rshape, fstat, &fcrit, &nodemean );

	for(max_spec=-1e20,i=0;i<nf;i++){		/* find max */
	  //spec[i] = sqrt(spec[i]);
	  if(spec[i] > max_spec)
	    max_spec = spec[i];
	}
	for(i=0;i<nf;i++)
	  // output is frequency power/max_power power  err
	  fprintf(stdout,"%20.15e %20.15e %20.15e %20.15e %20.15e\n",
		  freq[i],spec[i]/max_spec,spec[i],err[i],err[nf+i]);
	free(yk);free(wt);free(err);free(se);free(sk);free(fstat);free(spec);free(py);free(freq);
      }else{
	/* 
	   Lee and Park MT
	*/
	klen = jlpmtm_get_pow_2(nsample);
	klen = klen * pow((double) 2, (double) increase);
	nf = 1 + klen / 2;

	nyquist = 0.5 / dx;
	df = 2. * nyquist / klen;
	
	fprintf(stderr,"%s: MT LPMTM klen: %i nf: %i ny: %g df: %g\n",argv[0],klen,nf,nyquist,df);

	if((freq=(double *)calloc(klen,sizeof(double)))==NULL)
	  MEMERROR;
	if((spec=(double *)calloc(klen,sizeof(double)))==NULL)
	  MEMERROR;
	if((fstat=(double *)calloc(klen,sizeof(double)))==NULL)
	  MEMERROR;
	if((dof=(double *)calloc(klen,sizeof(double)))==NULL)
	  MEMERROR;
	nf = 1 + klen / 2;

	jlpmtm_do_mtap_spec(py,nsample, 1, kspec, hifac, 1, dx, spec, dof, fstat, klen);
	for(max_spec=-1e20,i=0;i<nf;i++){		/* find max */
	  spec[i]=spec[i]*spec[i];			/* square the output for consistencu */
	  if(spec[i] > max_spec)
	    max_spec = spec[i];
	}
	for(i=0;i < nf;i++){
	  freq[i] = df * i;
	  // output is frequency power/max_power power DOF Fvalues
	  fprintf(stdout,"%20.15e %20.15e %20.15e %20.15e %20.15e\n",
		  freq[i],spec[i]/max_spec,spec[i],dof[i],fstat[i]);
	}
	
	free(spec);free(dof);free(fstat);free(freq);
      }
    }
      
  }else{
    fprintf(stderr,"%s: too few data, %i\n",argv[0],n);
    exit(-1);
  }
}
COMP_PRECISION lin_interpolate ( COMP_PRECISION *x ,COMP_PRECISION *y, int n , COMP_PRECISION x1 )
{
  int i,j,n1;
  double fac,fac2;
  COMP_PRECISION tmp;
  n1 = n-1;
  j=0;
  while( (j < n1) && (x[j] < x1))
    j++;
  if(j == 0)
    j = 1;
  i = j - 1;
  
  fac = (x1 - x[i])/(x[j]-x[i]);
  fac2 = 1.0 - fac;
  tmp=  fac  * y[j] + fac2 * y[i];

  return (tmp);
}


/*

  
  from here only numerical recipes routines


*/

/*

  Lomb periodogram


*/
void period(COMP_PRECISION *x,COMP_PRECISION *y, int n,
	    COMP_PRECISION ofac, COMP_PRECISION hifac,
	    COMP_PRECISION *px,COMP_PRECISION *py,
	    int np,int *nout,int *jmax,
	    COMP_PRECISION *prob)
{
  int i,j;
  COMP_PRECISION ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh,
    sumsy,swtau,var,wtau,xave,xdif,xmax,xmin,yy;
  COMP_PRECISION arg,wtemp,*wi,*wpi,*wpr,*wr;
  
  wi=vector(1,n);
  wpi=vector(1,n);
  wpr=vector(1,n);
  wr=vector(1,n);
  *nout=0.5*ofac*hifac*n;
  if (*nout > np) nrerror("output arrays too short in period");
  avevar(y,n,&ave,&var);
  xmax=xmin=x[1];
  for (j=1;j<=n;j++) {
    if (x[j] > xmax) xmax=x[j];
    if (x[j] < xmin) xmin=x[j];
  }
  xdif=xmax-xmin;
  xave=0.5*(xmax+xmin);
  pymax=0.0;
  pnow=1.0/(xdif*ofac);
  for (j=1;j<=n;j++) {
    arg=TWOPID*((x[j]-xave)*pnow);
    wpr[j] = -2.0*SQUARE(sin(0.5*arg));
    wpi[j]=sin(arg);
    wr[j]=cos(arg);
    wi[j]=wpi[j];
  }
  for (i=1;i<=(*nout);i++) {
    px[i]=pnow;
    sumsh=sumc=0.0;
    for (j=1;j<=n;j++) {
      c=wr[j];
      s=wi[j];
      sumsh += s*c;
      sumc += (c-s)*(c+s);
    }
    wtau=0.5*atan2(2.0*sumsh,sumc);
    swtau=sin(wtau);
    cwtau=cos(wtau);
    sums=sumc=sumsy=sumcy=0.0;
    for (j=1;j<=n;j++) {
      s=wi[j];
      c=wr[j];
      ss=s*cwtau-c*swtau;
      cc=c*cwtau+s*swtau;
      sums += ss*ss;
      sumc += cc*cc;
      yy=y[j]-ave;
      sumsy += yy*ss;
      sumcy += yy*cc;
      wr[j]=((wtemp=wr[j])*wpr[j]-wi[j]*wpi[j])+wr[j];
      wi[j]=(wi[j]*wpr[j]+wtemp*wpi[j])+wi[j];
    }
    py[i]=0.5*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/var;
    if (py[i] >= pymax) pymax=py[(*jmax=i)];
    pnow += 1.0/(ofac*xdif);
  }
  expy=exp(-pymax);
  effm=2.0*(*nout)/ofac;
  *prob=effm*expy;
  if (*prob > 0.01) *prob=1.0-pow(1.0-expy,effm);
  free_vector(wr,1,n);
  free_vector(wpr,1,n);
  free_vector(wpi,1,n);
  free_vector(wi,1,n);
}

/*


  fast lomb periodogram
  
*/
void fasper(COMP_PRECISION *x,COMP_PRECISION *y,int n,
	    COMP_PRECISION ofac,COMP_PRECISION hifac,
	    COMP_PRECISION *wk1,COMP_PRECISION *wk2,
	    int nwk,int *nout,
	    int *jmax,COMP_PRECISION *prob)

{
  int j,k,ndim,nfreq,nfreqt;
  COMP_PRECISION ave,ck,ckk,cterm,cwt,den,df,effm,expy,fac,fndim,hc2wt;
  COMP_PRECISION hs2wt,hypo,pmax,sterm,swt,var,xdif,xmax,xmin;
  
  *nout=0.5*ofac*hifac*n;
  nfreqt=ofac*hifac*n*MACC;
  nfreq=64;
  while (nfreq < nfreqt) nfreq <<= 1;
  ndim=nfreq << 1;
  if (ndim > nwk) {
    fprintf(stderr,"workspaces too small in fasper: ndim: %i nwk: %i\n",
	    ndim,nwk);
    exit(-1);
  }
  avevar(y,n,&ave,&var);
  xmin=x[1];
  xmax=xmin;
  for (j=2;j<=n;j++) {
    if (x[j] < xmin) xmin=x[j];
    if (x[j] > xmax) xmax=x[j];
  }
  xdif=xmax-xmin;
  for (j=1;j<=ndim;j++) wk1[j]=wk2[j]=0.0;
  fac=ndim/(xdif*ofac);
  fndim=ndim;
  for (j=1;j<=n;j++) {
    ck=(x[j]-xmin)*fac;
    MOD(ck,fndim)
      ckk=2.0*(ck++);
    MOD(ckk,fndim)
      ++ckk;
    spread(y[j]-ave,wk1,ndim,ck,MACC);
    spread(1.0,wk2,ndim,ckk,MACC);
  }
  realft(wk1,ndim,1);
  realft(wk2,ndim,1);
  df=1.0/(xdif*ofac);
  pmax = -1.0;
  for (k=3,j=1;j<=(*nout);j++,k+=2) {
    hypo=sqrt(wk2[k]*wk2[k]+wk2[k+1]*wk2[k+1]);
    hc2wt=0.5*wk2[k]/hypo;
    hs2wt=0.5*wk2[k+1]/hypo;
    cwt=sqrt(0.5+hc2wt);
    swt=SIGN(sqrt(0.5-hc2wt),hs2wt);
    den=0.5*n+hc2wt*wk2[k]+hs2wt*wk2[k+1];
    cterm=SQUARE(cwt*wk1[k]+swt*wk1[k+1])/den;
    sterm=SQUARE(cwt*wk1[k+1]-swt*wk1[k])/(n-den);
    wk1[j]=j*df;
    wk2[j]=(cterm+sterm)/(2.0*var);
    if (wk2[j] > pmax) pmax=wk2[(*jmax=j)];
  }
  expy=exp(-pmax);
  effm=2.0*(*nout)/ofac;
  *prob=effm*expy;
  if (*prob > 0.01) *prob=1.0-pow(1.0-expy,effm);
}
/*


 */
void realft(COMP_PRECISION *data,int n,int isign)
{
  int i,i1,i2,i3,i4,np3;
  COMP_PRECISION c1=0.5,c2,h1r,h1i,h2r,h2i;
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
void four1(COMP_PRECISION *data,int nn,int isign)
{
  int n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  COMP_PRECISION tempr,tempi;
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
void spread(COMP_PRECISION y,COMP_PRECISION *yy,int n,COMP_PRECISION x,
	    int m)
{
  int ihi,ilo,ix,j,nden;
  static int nfac[11]={0,1,1,2,6,24,120,720,5040,40320,362880};
  COMP_PRECISION fac;
  if (m > 10) nrerror("factorial table too small in spread");
  ix=(int)x;
  if (x == (COMP_PRECISION)ix) yy[ix] += y;
  else {
    ilo=IMIN(IMAX((long)(x-0.5*m+1.0),1),n-m+1);
    ihi=ilo+m-1;
    nden=nfac[m];
    fac=x-ilo;
    for (j=ilo+1;j<=ihi;j++) fac *= (x-j);
    yy[ihi] += y*fac/(nden*(x-ihi));
    for (j=ihi-1;j>=ilo;j--) {
      nden=(nden/(j+1-ilo))*(j-ihi);
      yy[j] += y*fac/(nden*(x-j));
    }
  }
}
void nrerror(char *error_text)
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(-1);
}
COMP_PRECISION *vector(long nl,long nh)
{
  COMP_PRECISION *v;
  
  v=(COMP_PRECISION *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(COMP_PRECISION)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}
void free_vector(COMP_PRECISION *v,long nl,long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}
void avevar(COMP_PRECISION *data, long n,
	    COMP_PRECISION *ave,COMP_PRECISION *var)
{
  int j;
  COMP_PRECISION s,ep;
  
  for (*ave=0.0,j=1;j<=n;j++) *ave += data[j];
  *ave /= n;
  *var=ep=0.0;
  for (j=1;j<=n;j++) {
    s=data[j]-(*ave);
    ep += s;
    *var += s*s;
  }
  *var=(*var-ep*ep/n)/(n-1);
}

