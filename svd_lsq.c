/* 

sample program to solve a least squares problem by LAPACK SVD

this is written for clarity, not speed, and uses the high level driver
routines


$Id: svd_lsq.c,v 1.2 2011/01/02 22:01:45 becker Exp becker $



*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* 

lapack declarations

*/
/* needed for compiler */
#define sgelsd sgelsd_
#define sgelss sgelss_
#define dgelsd dgelsd_
#define dgelss dgelss_

/* SVD divide and conquer  */
extern void sgelsd(int *, int *, int *, float *,int *,
		   float *,int *,float *,float *,int *,
		   float *,int *, int *,int *);
extern void dgelsd(int *, int *, int *, double *,int *,
		   double *,int *,double *,double *,int *,
		   double *,int *, int *,int *);
/* regular SVD */
// SGELSS( M, N, NRHS, A, LDA, 
// B, LDB, S, RCOND, RANK, 
// WORK, LWORK, INFO )
extern void sgelss(int *, int *, int *, float *,int *,
		   float *,int *,float *,float *,int *,
		   float *,int *,int *);
extern void dgelss(int *, int *, int *, double *,int *,
		   double *,int *,double *,double *,int *,
		   double *,int *,int *);

/* 
   need to change here and above if we want double prec! 

*/
//#define DOUBLE_PRECISION

#ifdef DOUBLE_PRECISION
#define FFMT "%lf"
#define FFMT2 "%lf %lf"
#define CPREC double
#define LP_SVD_DC dgelsd
#define LP_SVD_REG dgelss
#else
#define CPREC float
#define FFMT "%f"
#define FFMT2 "%f %f"
#define LP_SVD_DC sgelsd
#define LP_SVD_REG sgelss
#endif

#define FEXIT {fprintf(stderr,"%s: file error\n",argv[0]);exit(-1);}
#define MEXIT {fprintf(stderr,"%s: memory error\n",argv[0]);exit(-1);}

int main(int argc, char **argv)
{
  int n,m,*iwork,rank,info;
  static int nrhs = 1;	    /*  */
  int mode = 1;		/* 1: divide and conquer 2: regular SVD */
  int binary = 0;	/* ASCII by default */
  int use_weights = 0;
  CPREC *a,*d,*work,*s,rcond,*acopy,dp,*dcopy,error,*w,tmp;
  FILE *in;
  int i,j,lwork,liwork,minmn,maxmn;
  char fname1[200],fname2[200],fname3[200],fname4[200],fname5[200];
  sprintf(fname1,"%s","a.dat");	/* matrix input */
  sprintf(fname2,"%s","b.dat");	/* data input */
  sprintf(fname3,"%s","x.dat");	/* solution out */
  sprintf(fname4,"%s","s.dat");	/* SVF out out */
  sprintf(fname5,"%s","bp.dat");	/* prediction output */
  if(argc < 4){
    fprintf(stderr,"%s m_data n_model rcond [%s] [%s] [%s] [%s] [mode, %i] [binary, %i] [use_weights, %i]\n",
	    argv[0],fname1,fname2,fname3,fname4,mode,binary,use_weights);
    fprintf(stderr,"\treads m by n matrix from %s,  m data from %s,\n",fname1,fname2);
    fprintf(stderr,"\tand solves the least squares problem |Ax-d| = min by SVD with rcond selection (see GELSD)\n");
    fprintf(stderr,"\tsolution vector will be in %s, SVD spectrum in %s\n",fname3,fname4);
    fprintf(stderr,"\tif rcond <=0, will use all SVs, else only if SV/SV_max >= rcond\n");
    fprintf(stderr,"\tmode=1: divide and conquer 2: regular SVD\n");
#ifdef DOUBLE_PRECISION
    fprintf(stderr,"\tbinary=0: ASCII input for all 1: A matrix is binary (double precision)\n");
#else
    fprintf(stderr,"\tbinary=0: ASCII input for all 1: A matrix is binary (single precision)\n");
#endif
    fprintf(stderr,"\tuse_weights: 0: read in data, 1: read and weight by data-weight pairs\n");
    exit(-1);
  }
  sscanf(argv[1],"%i",&m);
  sscanf(argv[2],"%i",&n);
  sscanf(argv[3],FFMT,&rcond);
  if(rcond <= 0)
    fprintf(stderr,"%s: using all singular values\n",argv[0]);
  if(rcond > 1){
   fprintf(stderr,"%s: rcond > 1 doesn't make sense\n",argv[0]);
   exit(-1);
  }

  if(argc>4)
    sprintf(fname1,"%s",argv[4]);
  if(argc>5)
    sprintf(fname2,"%s",argv[5]);
  if(argc>6)
    sprintf(fname3,"%s",argv[6]);
  if(argc>7)
    sprintf(fname4,"%s",argv[7]);
   if(argc>8)
    sscanf(argv[8],"%i",&mode);
  if(argc>9)
    sscanf(argv[9],"%i",&binary);
  if(argc>10)
    sscanf(argv[10],"%i",&use_weights);
  fprintf(stderr,"%s: allocating for m: %i by n: %i system, using rcond: %g binary: %i weights: %i\n",
	  argv[0],m,n,rcond,binary,use_weights);

  fprintf(stderr,"%s: reading from %s and %s, writing to %s and %s\n",
	  argv[0],fname1,fname2,fname3,fname4);
  
  /* allocate */
  minmn=(m<n)?(m):(n);
  maxmn=(m>n)?(m):(n);

  /*  */
  a=(CPREC *)malloc(sizeof(CPREC)*n*m);
  acopy=(CPREC *)malloc(sizeof(CPREC)*n*m);
  d=(CPREC *)malloc(sizeof(CPREC)*maxmn);
  dcopy=(CPREC *)malloc(sizeof(CPREC)*m);
  w=(CPREC *)malloc(sizeof(CPREC)*m);
  s=(CPREC *)malloc(sizeof(CPREC)*minmn);
  iwork=(int *)malloc(sizeof(int));
  if(mode == 1){		/* divide and conquer */
    work=(CPREC *)malloc(sizeof(CPREC));  
    /* memory inquiry */
    lwork=-1;
    LP_SVD_DC(&m,&n,&nrhs,a,&m,d,&maxmn,s,&rcond,&rank,work,&lwork,iwork,&info);

    liwork=iwork[0] + 500;
    lwork=(int)work[0];
    fprintf(stderr,"%s: divide and conquer: mem: iwork: %i lwork: %i\n",argv[0],liwork,lwork);

    /* properly allocate */
    work=(CPREC *)realloc(work,sizeof(CPREC)*lwork);  
    iwork=(int *)realloc(iwork,sizeof(int)*liwork);
  }else{
    /* regular */
    lwork = 3*minmn + ( 2*minmn > maxmn) ? (2*minmn):(maxmn);
    lwork *= 10;
    work=(CPREC *)malloc(sizeof(CPREC)*lwork);  

  }
  if(!a || !d || !s || !work || !acopy || !dcopy || !w)MEXIT;
  /* read in data */
  if(use_weights)
    fprintf(stderr,"%s: reading data and weights from %s\n",argv[0],fname2);
  else
    fprintf(stderr,"%s: reading data from %s\n",argv[0],fname2);
  in=fopen(fname2,"r");if(!in)FEXIT;
  if(use_weights){
    for(i=0;i < m;i++)
      if(fscanf(in,FFMT2,(d+i),(w+i))!=2)FEXIT; /* read in C style */
  }else{
    for(i=0;i < m;i++)
      if(fscanf(in,FFMT,(d+i))!=1)FEXIT; /* read in C style */
  }
	
  fclose(in);
  fprintf(stderr,"%s: read in data OK\n",argv[0]);
  memcpy(dcopy,d,sizeof(CPREC)*m); /* make a copy */
  
  /* read in a */
  if(binary){
#ifdef DOUBLE_PRECISION
    fprintf(stderr,"%s: reading matrix from %s (expecting double precision binary)\n",argv[0],fname1);
#else
    fprintf(stderr,"%s: reading matrix from %s (expecting single precision binary)\n",argv[0],fname1);
#endif
  }  else{
    fprintf(stderr,"%s: reading matrix from %s (expecting ASCII)\n",argv[0],fname1);
  }
  in = fopen(fname1,"r");if(!in)FEXIT;
  if(binary){
    for(i=0;i < m;i++)		/* m data rows */
      for(j=0;j < n;j++){		/* n parameter columns */
	/* Aij = a[j*m+i] */
	if(fread((a+j*m+i),sizeof(CPREC),1,in)!=1)
	  FEXIT;
      }
  }else{
    for(i=0;i < m;i++)		/* m data rows */
      for(j=0;j < n;j++){		/* n parameter columns */
	/* Aij = a[j*m+i] */
	if(fscanf(in,FFMT,(a+j*m+i))!=1)
	  FEXIT;
	/* 
	   read in FORTRAN style
	*/
      }
  }
  fclose(in);
  fprintf(stderr,"%s: read in matrix OK\n",argv[0]);
  /* make a copy */
  memcpy(acopy,a,sizeof(CPREC)*n*m);
  if(use_weights){		/* scale */
    for(i=0;i < m;i++){		/* m data rows */
      d[i] *= w[i];
      for(j=0;j < n;j++){		/* n parameter columns */
	a[j*m+i] *= w[i];
      }
    }
  }



  /* solve */
  //
  if(mode == 1){
    /* divide and conquor */
    LP_SVD_DC(&m,&n,&nrhs,a,&m,d,&maxmn,s,&rcond,&rank,work,&lwork,iwork,&info);
  }else{			/* regular SVD */
    LP_SVD_REG(&m,&n,&nrhs,a,&m,d,&maxmn,s,&rcond,&rank,work,&lwork,&info);
  }
  if(info != 0){
    fprintf(stderr,"%s: LAPACK SVD solver mode %i error code %i\n",
	    argv[0],mode,info);
    exit(-1);
  }
  fprintf(stderr,"%s: solver OK, effective rcond: %g yielded rank %i out of %i\n",
	  argv[0],rcond,rank,n);
			  
  fprintf(stderr,"%s: printing solution for X to %s\n",argv[0],fname3);
  in=fopen(fname3,"w");
  for(i=0;i < n;i++){
    //fprintf(stderr,"%12g\n",d[i]);
    fprintf(in,"%14.7e\n",d[i]);
  }
  fclose(in);

  /* singular values  */
  fprintf(stderr,"%s: printing SVD SVD/max to %s\n",argv[0],fname4);
  in=fopen(fname4,"w");
  for(i=0;i < n;i++){
    fprintf(in,"%14.7e %14.7e\n",s[i],s[i]/s[0]);
  }
  fclose(in);
  /* 
     aproximation of data B  
  */
  fprintf(stderr,"%s: printing b prediction to %s\n",argv[0],fname5);
  in=fopen(fname5,"w");
  error = 0;
  for(i=0;i<m;i++){
    dp = 0.0;
    for(j=0;j<n;j++)
      dp += acopy[j*m+i] * d[j];
    /* chi2 */
    tmp = (dp-dcopy[i])*w[i];	/* (y-yp)/sigma or multiplied by weights */
    error += tmp*tmp;
    fprintf(in,"%14.7e\n",dp);
  }
  fclose(in);
  fprintf(stderr,"%s: chi2 is %g, sqrt per parameter %g\n",
	  argv[0],error,sqrt(error)/(CPREC)n);
  
  free(a);free(d);free(work);free(s);free(iwork);free(acopy);free(dcopy);
  return 0;
}
