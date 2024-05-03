#include "two_layer.h"
/*
  based on twolayermodel.m in ShearWaveSplitting of the SplitLab package.
  
  following [Savage & Silver, 1994] 
  index 1 nominates the lower layer,
  index 2 represents the upper layer
  
  
  uses Kris Walker's implementation of Savage and Silver (1994) via
  Wuestefeld's Matlab implementation of SplitLab, converted to C
  
*/

/* pass model parameters fast azimuth fazi[2] in deg, dt[2] in deg,
   backzimuth in degrees returns effective fast axis, apparent delay
   time. 1 = bottom 2 = top
*/
void calc_two_layer(PREC *fazi, PREC *dt, PREC period, PREC bazi, 
		    PREC *faziapp, PREC *dtapp) /* output */
{
  PREC t,cos_t[2],sin_t[2],sum[2],phiapp;
  PREC ap,aps,Cc,Cs,alpha,sin_da,cos_da,da;
  PREC phi,a[2],sin_a[2],cos_a[2];
  PREC cos_alpha,sin_alpha,denom;
  int i;
  const PREC pi_fac = M_PI/180,pi2_fac=pi_fac*2;
  PREC omega_d2;
  omega_d2 = M_PI/period;	/* T = 8 second period , w = 2pi f, this is w/2  */
  /* special cases */
  if ((fazi[0] == fazi[1])||(fabs(fazi[0]-fazi[1]) == 90)||
      (fabs(fazi[0]-fazi[1]) == 180)||(fabs(fazi[0]-fazi[1]) ==270)||
      (fabs(fazi[0]-fazi[1]) == 360)){
    fazi[1] += EPS;   //THIS ENSURES APP IS NON-ZERO WITH PHI1=PHI2
    if (dt[0] == dt[1])
      dt[1] += EPS;
  }
  //calculate expressions
  for(i=0;i<2;i++){
    phi = fazi[i] - bazi;	/* fast azimuth relative to
				   back-azimuth */
    a[i] = pi2_fac*phi;
    SCFUNC(a[i],(sin_a+i),(cos_a+i));
    /* general model parameters */
    t = omega_d2*dt[i];
    SCFUNC(t,(sin_t+i),(cos_t+i));
  }
  da = a[1]-a[0];
  SCFUNC(da,&sin_da,&cos_da);

  Cc  = cos_t[0]*sin_t[1]*cos_a[1] + cos_t[1]*sin_t[0]*cos_a[0];
  Cs  = cos_t[0]*sin_t[1]*sin_a[1] + cos_t[1]*sin_t[0]*sin_a[0];

  ap  =  cos_t[0]*cos_t[1] - sin_t[0]*sin_t[1]*cos_da;
  aps = -sin_t[0]*sin_t[1] * sin_da;

  /*  */
  denom = aps*ap+Cs*Cc;
  if(denom == 0)
    denom = EPS;		/* should we do this differently,
				   using atan2? */
  alpha = atan((aps*aps + Cs*Cs)/denom);
  SCFUNC(alpha,&sin_alpha,&cos_alpha);

  /*  */
  phiapp = alpha/pi2_fac;					/* phi =  alpha/2 */
  *faziapp = phiapp + bazi;				/* apparent fast azimuth */
  if(*faziapp>180)
    *faziapp -= 180;
  
  sum[0] = ap*sin_alpha - aps * cos_alpha;
  sum[1] = Cs*cos_alpha - Cc*sin_alpha;
  if(fabs(sum[0]) > EPS)
    *dtapp = atan(Cs/sum[0])/omega_d2;
  else
    *dtapp = atan(aps/sum[1])/omega_d2;
  if(*dtapp < 0){
    *dtapp = -(*dtapp);
    *faziapp += 90;
  }
  if(*faziapp > 180)		/* return in -90...90 range */
    *faziapp -= 180;
  if(*faziapp > 90)
    *faziapp -=180;
}

/* 
   convert short storage solution to full solution 
*/
void sli2sl(struct sli sol_short,struct sl *sol)
{
  int i;
  sol->cost = sol_short.cost;
  for(i=0;i<2;i++){
    sol->fazi[i] = -90 + TL_DA * (((PREC)sol_short.iazi[i])+0.5);
    sol->dt[i]   =       TL_DT * (((PREC)sol_short.idt[i])+0.5);
  }
}
