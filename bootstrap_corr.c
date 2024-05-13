#include "correl_nr.h"
#define COMP_PRECISION double

COMP_PRECISION corr(COMP_PRECISION *,COMP_PRECISION *, int );
void replicate(COMP_PRECISION *, COMP_PRECISION *, COMP_PRECISION *, COMP_PRECISION *, int , long  *,int, int *); /*  */
COMP_PRECISION ran3(long *);
int random_int(long *, int );

void calc_avg(COMP_PRECISION *dx,COMP_PRECISION *,int ,COMP_PRECISION *, COMP_PRECISION *,int,long int *);

/* 
   
   compute bootstrap errors on a correlation

 */
#define NBOOT 100000

int main(int argc, char **argv)
{
  int n;
  long int seed = -1;
  COMP_PRECISION *x, *y;
  COMP_PRECISION mean, std;
  
  /* init */
  ran3(&seed);
  /* allocate */

  x=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
  y=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
  /* read in */
  n=0;
  while(fscanf(stdin,"%lf %lf",(x+n),(y+n))==2){
    n++;
    x=(COMP_PRECISION *)realloc(x,(n+1)*sizeof(COMP_PRECISION));
    y=(COMP_PRECISION *)realloc(y,(n+1)*sizeof(COMP_PRECISION));
  }
#ifdef USE_PEARSON
  fprintf(stderr,"%s: read %i samples for Pearson correlation\n",argv[0],n);
#else
  fprintf(stderr,"%s: read %i samples for Spearman rank correlation\n",argv[0],n);
#endif
  //calc_avg(x,y,n,&mean,&std,1,&seed);
  //  fprintf(stderr,"%s: jacknife:                     %11g %11g\n",argv[0],mean,std);

  calc_avg(x,y,n,&mean,&std,2,&seed);
  fprintf(stderr,"%s: bootstrap (%8i samples): %11g %11g\n",argv[0],NBOOT,mean,std);
  printf("%24.7e %14.7e\n",mean,std);
  
  
  exit(-1);
}
/* 
   1: jacknife (usually, bad)
   2: bootstrap (usually, good)

 */
void calc_avg(COMP_PRECISION *x,COMP_PRECISION *y,int n,COMP_PRECISION *gmean,COMP_PRECISION *std,int mode,long int *seed)
{
  double sum,mean;
  int nused,nrep,i,nnew;
  COMP_PRECISION *rx,*ry;
  COMP_PRECISION d,zd,probd,probrs,rtmp;
  double *r;
  /* for replication */
  rx=(COMP_PRECISION *)malloc(n*sizeof(COMP_PRECISION));
  ry=(COMP_PRECISION *)malloc(n*sizeof(COMP_PRECISION));
  
  if(mode==1){			/* jacknife */
    nrep = n;
  }else if(mode==2){			/* bootstrap */
    nrep = NBOOT;
  }else{
    fprintf(stderr,"mode %i undefined\n",mode);
    exit(-1);
  }
  r=(double *)malloc(nrep*sizeof(double));

  sum =0;
  nused=0;
  for(i=0;i < nrep;i++){		/* sample */
    if(mode==1){			/* jacknife */
      replicate(x,y,rx,ry,n,seed,-(1+i),&nnew);
    }else{			/* bootstrap */
      replicate(x,y,rx,ry,n,seed,1,&nnew);
    }
       
#ifdef USE_PEARSON
    rtmp = corr(rx,ry,nnew);
#else
    spear ((rx-1),(ry-1),(unsigned long)nnew,&d,&zd, &probd, &rtmp,&probrs);
#endif
    if(finite(rtmp)){
      sum += (double)rtmp;
      r[nused] = (double)rtmp;
      nused++;
    }
  }
  mean = sum / (double)nused;
  *gmean = (COMP_PRECISION)mean;
  
  sum = 0;
  for(i=0;i < nused;i++){
    r[i] -= mean;
    sum += r[i]*r[i];
  }


  *std = (COMP_PRECISION)sqrt(sum/((double)(nused-1)));

  free(rx);free(ry);free(r);
}

/* 
   1: provide a random sampling with replicates 
   2: take out one random sample
   
*/
void replicate(COMP_PRECISION *x, COMP_PRECISION *y, COMP_PRECISION *rx, COMP_PRECISION *ry, int n, long *seed,
	       int mode,int *nnew)
{
  int i,j,omit;
  if(mode == 1){/* resample all, allowing replicates */
    for(i=0;i < n;i++){
      j = random_int(seed,n);
      //fprintf(stderr,"%i\n",j);
      rx[i] = x[j];
      ry[i] = y[j];
    }
    *nnew = n;
  }else if(mode == 2){		/* take out one random sample */
    omit = random_int(seed,n);
    for(i=j=0;i < n;i++){
      if(i != omit){
	rx[j] = x[i];
	ry[j] = y[i];
	j++;
      }
    }
    *nnew = j;
  }else if(mode < 0){		/* take out sample -mode-1 */
    omit = -mode-1;
    if((omit<0)||(omit>n-1)){
      fprintf(stderr,"index %i error mode %i\n",omit,mode);
      exit(-1);
    }
    for(i=j=0;i < n;i++){
      if(i != omit){
	rx[j] = x[i];
	ry[j] = y[i];
	j++;
      }
    }
    *nnew = j;
  }else{
    fprintf(stderr,"error mode %i\n",mode);
    exit(-1);
  }
  //fprintf(stderr,"\n");
}

COMP_PRECISION corr(COMP_PRECISION *x,COMP_PRECISION *y, int n)
{
  COMP_PRECISION mx,my,s1,s2,s3,tmp,dx,dy;
  int i;
  mx=0.0;
  my=0.0;
  for(i=0;i<n;i++){
    mx += x[i];
    my += y[i];
  }
  mx /= (COMP_PRECISION)n;
  my /= (COMP_PRECISION)n;
  s1=s2=s3=0.0;
  for(i=0;i<n;i++){
    dx = x[i] - mx;
    dy = y[i] - my;
    s1 += (dx*dy);
    s2 += (dx*dx);
    s3 += (dy*dy);
  }
  tmp = sqrt(s2) * sqrt(s3);
  return s1/tmp;

}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

COMP_PRECISION ran3(long *idum)
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  
  if (*idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

/* return integer between 0....n-1 */
int random_int(long  *seed, int n)
{
  int j;
  j=0;
  j = (int)(ran3(seed)*(COMP_PRECISION)n);
  //fprintf(stderr,"%i %g\n",j,(COMP_PRECISION)j/(COMP_PRECISION)(n-1));
  if(j==n)
    j--;
  return j;
}
