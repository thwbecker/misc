#include "fit_exp.h"
/* 


 */

#define NMAX 4
#define MODE_MAX 2
int main(int argc, char **argv)
{
  double par_init[MODE_MAX][NMAX]={{1,1000,0,0},{0,0,1,1000}},par[NMAX]; 
  int m = 0,n,mode,i;
  double *t ;
  double *y ;
  double (*fit_func)();
  double norm;
  const int verbose = 0;
  t=(double *)malloc(sizeof(double));
  y=(double *)malloc(sizeof(double));

  //in=fopen("tmp.dat","r");
  while(fscanf(stdin,"%lf %lf\n",(t+m),(y+m))==2){
    m++;
    t=(double *)realloc(t,sizeof(double)*(m+1));
    y=(double *)realloc(y,sizeof(double)*(m+1));
  }
  fprintf(stderr,"%s: read %i data points\n",argv[0],m);
  //fclose(in);
  for(mode=0;mode < MODE_MAX;mode++){
    for(i=0;i<NMAX;i++)
      par[i] = par_init[mode][i];
	  
    switch(mode){
    case 0:
      n=2;fit_func=fit_exp1;
      break;
    case 1:
      n=4;fit_func=fit_exp2;
      break;
    }
    norm = fit_exp_driver(m,t,y,par,verbose,3,fit_func);
    printf("%i %g ",mode,norm);
    for(i=0;i<n;i++)
      printf("%g ",par[i]);
    printf("\n");
  }
  free(t);
  free(y);

  return 0;
}


double fit_exp1( double t, const double *p )
{
  //return p[0] + p[1] * (1-exp(-t/p[2]));
  return p[0] * (1-exp(-t/(p[1]*p[1])));  
}
double fit_exp2( double t, const double *p )
{
  return p[0] + p[1] * t + p[2]*(1-exp(-t/(p[3]*p[3])));  
  //return p[0] + p[1] * (1- p[3]* exp(-t/p[2]) + (1-p[3])*exp(-t/p[4]));
}



double fit_exp_driver(int m, double *t, double *y, double *par, int verbose,
		    int n,double(*fit_func)(double,const double *))
{
  lm_control_struct control = lm_control_double;
  lm_status_struct status;
  int i;
  control.verbosity = 0;
  if(verbose)
    printf( "Fitting ...\n" );

  /* now the call to lmfit */
  lmcurve( n, par, m, t, y, fit_func, &control, &status);
  
  if(verbose){
    printf( "Results:\n" );
    printf( "status after %d function evaluations:\n  %s\n",
	    status.nfev, lm_infmsg[status.outcome] );
  
    printf("obtained parameters:\n");
    for ( i = 0; i < n; ++i)
      printf("  par[%i] = %12g\n", i, par[i]);
    printf("obtained norm:\n  %12g\n", status.fnorm );
    
    printf("fitting data as follows:\n");
    for ( i = 0; i < m; ++i)
      printf( "  t[%2d]=%4g y=%6g fit=%10g residue=%12g\n",
	      i, t[i], y[i], 
	      fit_func(t[i],par), y[i] - fit_func(t[i],par) );
  }
  return status.fnorm;
}




