/* 
   
   compute mean and cov from event times

 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NR 50000		/* bootstrap */



void calm_and_s(double *,double *,double *,double *,int );
void calm_and_s_from_vec(float *,double *,double *,int );

int main(int argc, char **argv)
{
  double *e,tsm,tss,cov;
  double ts_mean,ts_std,cov_mean,cov_std;
  float *vt,*vc;
  double *dt,*dtc;
  int i,k,n,j,nel,ind;
  int rmode = 0;		/* 0: bootstrap 1: jacknife */
  if(argc > 1)
    sscanf(argv[1],"%i",&rmode);
  
  e=(double *)malloc(sizeof(double));
  n=0;
  while(fscanf(stdin,"%lf",(e+n))==1){
    n++;
    e=(double *)realloc(e,sizeof(double)*(n+1));
  }
  fprintf(stderr,"%s: read %i event times, error mode: %i\n",argv[0],n,rmode);
  dt = (double *) malloc(sizeof(double)*(n-1));
  dtc = (double *) malloc(sizeof(double)*(n-1));

  if((!dt)||(!dtc)){
    fprintf(stderr,"%s: mem error %i\n",argv[0],n-1);
    exit(-1);
  }
  for(i=1;i<n;i++)
    dt[i-1] = e[i]-e[i-1];
  n--;
  free(e);

  calm_and_s(dt,&tsm,&tss,&cov,n);
  fprintf(stderr,"%s: Ts: %7.1f cov: %.2f N_all: %i\n",argv[0],tsm,cov,n+1);
  if(rmode==0)
    nel = NR;
  else
    nel = n;
  vt=(float *)malloc(sizeof(float)*nel);
  vc=(float *)malloc(sizeof(float)*nel);
  if((!vt)||(!vc)){
    fprintf(stderr,"%s: mem error %i\n",argv[0],nel);
    exit(-1);
  }
  if(rmode == 0){		/* bootstrap */
    
    for(i=0;i<NR;i++){
      for(j=0;j<n;j++){
	ind = (int)((float)random()/(float)RAND_MAX*(float)n);
	if(ind==n)
	  ind--;
	dtc[j] = ind;
      }
      calm_and_s(dtc,&tsm,&tss,&cov,n);
      //fprintf(stderr,"%s: Ts: %7.1f cov: %.2f N_del: %i\n",argv[0],tsm,cov,n+1-ndel);
      vt[i] = tsm;vc[i] = cov;	/* store recurrence time and cov of this sample */
    }
  }else if(rmode == 1){		/* jacknife */
    for(i=0;i < nel;i++){
      for(k=j=0;j < n;j++){
	if(j != i){
	  dtc[k] = dt[j];
	  k++;
	}
      }
      calm_and_s(dtc,&tsm,&tss,&cov,n-1);
      vt[i] = tsm;vc[i] = cov;
    }
  }else{
    fprintf(stderr,"%s: error mode %i undefined\n",argv[0],rmode);
    exit(-1);
  }
    
  calm_and_s_from_vec(vt,&ts_mean,&ts_std,nel);
  calm_and_s_from_vec(vc,&cov_mean,&cov_std,nel);
  fprintf(stderr,"%s: mean Ts: %7.1f %.2f cov: %.2f %.3f based on %i resamples\n",
	  argv[0],ts_mean,ts_std,cov_mean,cov_std,nel);
  fprintf(stdout,"%20.7e %20.7e %20.7e %20.7e\n",ts_mean,ts_std,cov_mean,cov_std);
  free(vt);free(vc);free(dt);free(dtc);
  return 0;
}
/* calc mean and std of intervals */
void calm_and_s(double *dt,double *tsm,double *tss,double *cov,int n)
{
  int i;
  double tmp;
  if(n<2){
    fprintf(stderr,"calm_and_s: n too small %i\n",n);
    exit(-1);
  }
  *tsm = *tss = 0;
  for(i=0;i<n;i++)
    *tsm += dt[i];
  *tsm /= (double)n;
  for(i=0;i<n;i++){
    tmp = dt[i] - (*tsm);
    *tss += tmp*tmp;
  }
  *tss = sqrt(*tss/(n-1));
  *cov = (*tss)/(*tsm);
}
/* only mean and std on vector */
void calm_and_s_from_vec(float *x,double *tsm,double *tss,int n)
{
  int i;
  double tmp;
  if(n<2){
    fprintf(stderr,"calm_and_s_from_vec: n too small %i\n",n);
    exit(-1);
  }
  *tsm = *tss = 0;
  for(i=0;i<n;i++)
    *tsm += x[i];
  *tsm /= (double)n;
  for(i=0;i<n;i++){
    tmp = x[i] - (*tsm);
    *tss += tmp*tmp;
  }
  *tss = sqrt(*tss/(n-1));
}

