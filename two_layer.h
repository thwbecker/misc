#define _GNU_SOURCE    
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>

void sincos(double, double *, double *);
void sincosf(float, float *, float *);


#define SINGLE_PREC

#ifdef SINGLE_PREC
#define PREC float		/* single prec */
#define SCFUNC sincosf		/* sincos function */
#define FFMT "%f"
#define F5FMT "%f %f %f %f %f"
#define EPS 1e-7
#else  /* double */
#define PREC double
#define SCFUNC sincos		/* sincos function */
#define FFMT "%lf"
#define F5FMT "%lf %lf %lf %lf %lf"
#define EPS 1e-15

#endif

/* global discretization */
//#define TL_DA 5.0		/* 8.32 is the mean error */
//#define TL_DT 0.1		/* 0.25 is the mean error */

#define TL_DA 2.5		/* 8.32 is the mean error */
#define TL_DT 0.05		/* 0.25 is the mean error */

/* split structure */
struct sp{
  PREC phi,wphi,dt,wdt,ba;
};
/* misfit structure */
struct mf{
  PREC da;
  PREC dt;
};
/* solution structure */
struct sl{
  PREC fazi[2];
  PREC dt[2];
  PREC cost;
};
/* compressed solution structure*/
struct sli{
  unsigned char iazi[2];
  unsigned char idt[2];
  PREC cost;
};
struct nb{
  PREC weight;
  unsigned int id;
};
/* station structure */
struct st{
  char name[300];
  PREC lon, lat,*dist;
  int nsol,nnb;
  struct nb *neighbor;
  struct sli *sol;
  struct sl sol2l,sol1l;
  int nobs;
  PREC wfazi,period;
};
/* model structure */
struct sts{
  int nstat;
  struct st *station;

  PREC aw,lw,dw,tlw;
  long rseed;
};
/* individual solution structure */
struct ind{
  int *gene;
  PREC fitness,cost,scost,dcost;
};
void calc_fitness(struct ind  *, struct sts);
void calc_two_layer(PREC *,PREC *,PREC,PREC , PREC *, PREC *); 
PREC cdirdiff(PREC ,PREC, int );
PREC cost_function(struct mf,PREC,int);
void add_sol(struct sl *, struct sl ,int *,int);
int compare_soli(struct sli *, struct sli *);
int compare_ind(struct ind *, struct ind *);
void calc_chi2(struct mf *,struct sp *, int ,struct mf,struct sl ,PREC );
PREC median(PREC *, int );
PREC nr_select(int ,int ,PREC *);
void sli2sl(struct sli ,struct sl *);
PREC distance(PREC ,PREC ,PREC ,PREC );
PREC distance_rad(PREC,PREC ,PREC ,PREC , PREC, PREC );
PREC calc_cost(struct sts , int *,int, PREC *, PREC *);
void perturb_sol_gene(struct sts *, int *, int *,PREC);
PREC ran2(long *);
PREC ran3(long *);
void mate(struct ind , struct ind ,struct ind *, struct sts *,PREC );
