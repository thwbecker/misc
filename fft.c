#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/*

  rudimentary FFT routine, see also period

  

  $Id: fft.c,v 1.1 2005/02/23 03:04:21 becker Exp becker $



 */
#ifndef SQUARE	
static double sqrargsdf;
#define SQUARE(x) (((sqrargsdf=(x))== 0.0) ? (0.0) : (sqrargsdf*sqrargsdf))
#endif
#define ABS(x) (((x)>=0)?((x)):(-(x)))

#define DEF_MODE 1
#define COMP_PRECISION double

void phelp(char *);
COMP_PRECISION *vector(long ,long );
void free_vector(COMP_PRECISION *,long ,long );
void four1(COMP_PRECISION *,unsigned long ,int );
void spctrm(FILE *,COMP_PRECISION *,int ,int , int );
int nextpwrtwo(int );

#define WINDOW(j, hn) (1.0) // general window for abs values
/* Bartlett window for periodogram */
#define SPEC_WINDOW(j,a,b) (1.0-fabs((((j)-1)-(a))*(b)))       

int main(int argc, char *argv[] )
{
  int i,j,m,hn,n,mode;
  COMP_PRECISION *f,ftmp,delta;
  FILE *in;
  char filename[200];

  delta=1.0;
  mode=DEF_MODE;

  if(argc > 1 && strcmp(argv[1],"-h") == 0 )argc=-9;
  switch(argc){
  case 2:{
    sprintf(filename,"%s",argv[1]);
    break;
  }
  case 3:{
    sprintf(filename,"%s",argv[1]);
    sscanf(argv[2],"%lf",&delta);
    break;
  }
  case 4:{
    sprintf(filename,"%s",argv[1]);
    sscanf(argv[2],"%lf",&delta);
    sscanf(argv[3],"%i",&mode);
    break;
  }
  default:{
    phelp(argv[0]);
    break;
  }}
  
  if((in=fopen(argv[1],"r"))==NULL){
    fprintf(stderr,"%s: couldn't open %s\n",argv[0],argv[1]);
    exit(-1);
  }
  /* count data */
  m=0;
  while(fscanf(in,"%lf\n",&ftmp)==1)
    m++;
  rewind(in);
  fprintf(stderr,"%s: number of data points: %i, Delta T: %g\n",argv[0],m,delta);
  
  
  
  switch(mode){
  case 1:{
    // periodogram
    fprintf(stderr,"%s: periodogram output\n",argv[0]);
    n=nextpwrtwo(m);// check if power of 2
    if(n != m){
      fprintf(stderr,"need power of two data points for spectral estimate\n");
      exit(-1);
    }
    hn= n / 2;
    if((f=(COMP_PRECISION *)calloc(n,sizeof(COMP_PRECISION)))==NULL){
      fprintf(stderr,"memory allocation error\n");
      exit(-1);
    }
    spctrm(in,f-1,hn/2,1,0);
    for(i=0;i< hn;i++)
      printf("%g %20.15e\n",
	     (COMP_PRECISION)i/((COMP_PRECISION)hn*delta*(COMP_PRECISION)n),
	     f[i]);
    break;
  }
  case 0:{
    //
    // simple fft output
    // 
    fprintf(stderr,"%s: fft output\n",argv[0]);
    n=nextpwrtwo(m);// check if power of 2
    if(n-m)
      fprintf(stderr,"%s: padding %i zeros to get 2^%i=%i\n",argv[0],n-m,i,n);
    hn= n / 2;
    if((f=(COMP_PRECISION *)calloc(2*n,sizeof(COMP_PRECISION)))==NULL){
      fprintf(stderr,"memory allocation error\n");
      exit(-1);
    }
    for(i=0;i<2*m;i+=2)
      fscanf(in,"%lf",(f+i));
    four1(f-1,n,1);
    // output 
    for(i=0;i<2*n;i+=2)
      printf("%20.15e\n%20.15e\n",f[i],f[i+1]);
    break;
  }
  //
  // absolute power of FFT output
  //
  case 2:{
    fprintf(stderr,"%s: absolute value of fft output\n",argv[0]);
    n=nextpwrtwo(m);// check if power of 2
    if(n-m)
      fprintf(stderr,"%s: padding %i zeros to get 2^%i=%i\n",argv[0],n-m,i,n);
    hn= n / 2;
    if((f=(COMP_PRECISION *)calloc(2*n,sizeof(COMP_PRECISION)))==NULL){
      fprintf(stderr,"memory allocation error\n");
      exit(-1);
    }
    if(WINDOW(0,hn)!=1.0){
      fprintf(stderr,"%s using taper\n",argv[0]);
    }
    for(j=i=0;i<m*2;i += 2,j++){
      fscanf(in,"%lf",(f+i));
      f[i] *= WINDOW(j,hn);
    }
    four1(f-1,n,1);

    ftmp=1.0/(((COMP_PRECISION)n)*delta);

    for(i=0;i<n+1;i += 2)
      printf("%20.15e %20.15e\n",
	     (COMP_PRECISION)i*ftmp,
	     sqrt(SQUARE(f[i])+SQUARE(f[i+1])));
    break;
  }
  default:{
    fprintf(stderr,"%s: mode %i is undefined\n",argv[0],mode);
    exit(-1);
    break;
  }}
  
  fclose(in);
  return 0;
}


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(COMP_PRECISION *data,unsigned long nn,int isign)
{
  unsigned long n,mmax,m,j,istep,i;
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
#undef SWAP
#undef WINDOW


void spctrm(FILE *fp,COMP_PRECISION *p,
	    int m,int k, int ovrlap)
{
  int mm,m44,m43,m4,kk,joffn,joff,j2,j;
  COMP_PRECISION w,facp,facm,*w1,*w2,
    sumw=0.0,den=0.0;

  mm=m+m;
  m43=(m4=mm+mm)+3;
  m44=m43+1;
  w1=vector(1,m4);
  w2=vector(1,m);
  facm=m;
  facp=1.0/m;
  for (j=1;j<=mm;j++) sumw += SQUARE(SPEC_WINDOW(j,facm,facp));
  for (j=1;j<=m;j++) p[j]=0.0;
  if (ovrlap)
    for (j=1;j<=m;j++) fscanf(fp,"%lf",&w2[j]);
  for (kk=1;kk<=k;kk++) {
    for (joff = -1;joff<=0;joff++) {
      if (ovrlap) {
	for (j=1;j<=m;j++) w1[joff+j+j]=w2[j];
	for (j=1;j<=m;j++) fscanf(fp,"%lf",&w2[j]);
	joffn=joff+mm;
	for (j=1;j<=m;j++) w1[joffn+j+j]=w2[j];
      } else {
	for (j=joff+2;j<=m4;j+=2)
	  fscanf(fp,"%lf",&w1[j]);
      }
    }
    for (j=1;j<=mm;j++) {
      j2=j+j;
      w=SPEC_WINDOW(j,facm,facp);
      w1[j2] *= w;
      w1[j2-1] *= w;
    }
    four1(w1,mm,1);
    p[1] += (SQUARE(w1[1])+SQUARE(w1[2]));
    for (j=2;j<=m;j++) {
      j2=j+j;
      p[j] += (SQUARE(w1[j2])+SQUARE(w1[j2-1])
	       +SQUARE(w1[m44-j2])+SQUARE(w1[m43-j2]));
    }
    den += sumw;
  }
  den *= m4;
  for (j=1;j<=m;j++) p[j] /= den;
  free_vector(w2,1,m);
  free_vector(w1,1,m4);
}



#define NR_END 1
#define FREE_ARG char*


COMP_PRECISION *vector(long nl,long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  COMP_PRECISION *v;
  v=(COMP_PRECISION *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(COMP_PRECISION)));
  if (!v){fprintf(stderr,"allocation failure in vector");exit(-1);}
  return v-nl+NR_END;
}

void free_vector(COMP_PRECISION *v,long nl,long nh)
     /* free a float vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}
void phelp(char *program)
{
  fprintf(stderr,"%s filename [Delta T, 1.0] [mode, %i]\n",program,DEF_MODE);
  fprintf(stderr,"\tmode 0: FFT output\n");
  fprintf(stderr,"\tmode 1: prediodogram output\n");
  fprintf(stderr,"\tmode 2: absolute value of FFT output\n");
  exit(-1);

}
int nextpwrtwo(int i)
{
  COMP_PRECISION y;
  y = ceil(log((COMP_PRECISION)i)*1.44269504088896341);
  return (int)pow(2,y);
}

