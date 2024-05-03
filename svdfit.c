/* 
   do SVD fit using plain numerical recipes for testing
   purposes. modified in that the svdfit routine takes the UNSCALED
   matrix A instead of a function afunc.  A is afunc evaluated at the
   X locations

   $Id: svdfit.c,v 1.2 2005/04/16 01:30:12 becker Exp becker $
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nr_defines.h"
void svdfit(COMP_PRECISION *,COMP_PRECISION *,int ,COMP_PRECISION *,int ,
	    COMP_PRECISION **,COMP_PRECISION **,COMP_PRECISION *,COMP_PRECISION *,
	    COMP_PRECISION **,COMP_PRECISION);

#define READERR {fprintf(stderr,"%s: read error\n",argv[0]);exit(-1);}
#define TOL_DEF 1.0e-7

int main(int argc, char **argv)
{
  COMP_PRECISION *y,*sig,*a,**u,**v,chisq,**cvm,*w,**amat,thres=-1;
  int i,j,ndata,m,fmode,use_sigma=0;
  char afile[300],yfile[300];
  FILE *in;
  if(argc < 3){
    fprintf(stderr,"%s ndata mmodel [a.dat] [y.dat] [y_has_uncertainty, 0] [thres, -1]\n",argv[0]);
    
    exit(-1);
  }else{
    sscanf(argv[1],"%i",&ndata);
    sscanf(argv[2],"%i",&m);
    if(argc>3)
      sprintf(afile,argv[3]);
    else
      sprintf(afile,"a.dat");
    if(argc>4)
      sprintf(yfile,argv[4]);
    else
      sprintf(yfile,"y.dat");
    if(argc > 5)
      sscanf(argv[5],"%i",&use_sigma);
    if(argc > 6)
      sscanf(argv[6],"%lf",&thres);
  }
  
  y=nr_vector(1,ndata);sig=nr_vector(1,ndata);a=nr_vector(1,m);
  w=nr_vector(1,m);cvm=nr_matrix(1,m,1,m);
  amat=nr_matrix(1,ndata,1,m);
  u=nr_matrix(1,ndata,1,m);
  v=nr_matrix(1,m,1,m);

  in = fopen(yfile,"r");
  if(!in){
    fprintf(stderr,"%s: error: %s not found\n",
	    argv[0],yfile);
    exit(-1);
  }
  if(use_sigma){
    /* we have uncertaintites */
    for(i=1;i<=ndata;i++){
      if(fscanf(in, TWO_FLT_FORMAT,(y+i),(sig+i))!=2)
	READERR;
    }
    fprintf(stderr,"%s: %i data with uncertainties ok\n",
	    argv[0],ndata);
  }else{
    /* no uncertainties */
    for(i=1;i<=ndata;i++){
      if(fscanf(in,FLT_FORMAT,(y+i))!=1)
	READERR;
      sig[i] = 1.0;
    }
    fprintf(stderr,"%s: %i data without uncertainties ok\n",
	    argv[0],ndata);
  }
  fclose(in);
  /* 
     read in amat 
  */
  in=fopen(afile,"r");
  if(!in){
    fprintf(stderr,"%s: couldn't find %s\n",
	    argv[0],afile);
    exit(-1);
  }
  for(i=1;i<=ndata;i++){
    for(j=1;j<=m;j++){
      if(fscanf(in,FLT_FORMAT,&amat[i][j])!=1)
	READERR;
    }
  }
  fclose(in);
  fprintf(stderr,"%s: read %i data and %i parameters ok\n",
	  argv[0],ndata,m);
  /* 
     solve 
  */
  svdfit(y,sig,ndata,a,m,u,v,w,&chisq,amat,thres);
  fprintf(stderr,"%s: chi2: %g rchi2: %g\n",argv[0],chisq, chisq/((ndata-m>0)?(ndata-m):(1)));
  /* calculate uncertainty */
  nr_svdvar(v,m,w,cvm);
  /* output of solution parameters and uncertainties */
  for(i=1;i<=m;i++)
    printf("%20.7e %20.7e\n",a[i],sqrt(cvm[i][i]));

  return 0;
}
/* amat[i][j] is the test function  evaluated at x[i] for the 
   j-th parameter 

   thres: if -1, use all singular values, else give fraction of largest

*/

void svdfit(COMP_PRECISION *y,COMP_PRECISION *sig,int ndata,COMP_PRECISION *a,int ma,
	    COMP_PRECISION **u,COMP_PRECISION **v,COMP_PRECISION *w,COMP_PRECISION *chisq,
	    COMP_PRECISION **amat,COMP_PRECISION thres)
{
  int j,i,iz;
  COMP_PRECISION wmax,tmp,sum,*b;
  if(thres < 0)
    thres = TOL_DEF;
  b=nr_vector(1,ndata);

  for (i=1;i<=ndata;i++) {
    tmp=1.0/sig[i];
    for (j=1;j<=ma;j++) 
      u[i][j] = amat[i][j]*tmp;
    b[i]=y[i]*tmp;
  }
  /* do decomposition */
  nr_svdcmp(u,ndata,ma,w,v);
  /* find max */
  wmax=0.0;
  for (j=1;j<=ma;j++)
    if (w[j] > wmax) wmax=w[j];

  fprintf(stderr,"svdfit: thres: %g ",thres);
  thres *= wmax;

  /* truncate */
  iz=0;
  for (j=1;j<=ma;j++){
    if (w[j] < thres) {
      w[j]=0.0;
      iz++;
    }
  }
  fprintf(stderr,"lead to %i singular values out of %i be set to zero\n",
	  iz,ma);
  /* backsub */
  nr_svbksb(u,w,v,ndata,ma,b,a);
  
  *chisq=0.0;
  for (i=1;i<=ndata;i++) {
    for (sum=0.0,j=1;j<=ma;j++) 
      sum += a[j]*amat[i][j];
    *chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
  }
  nr_free_vector(b,1,ndata);
}
#undef TOL_DEF
