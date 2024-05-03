#include "determine_netr.h"

void determine_netr(float *lon,float *lat,
		    float *velt,float *velp,
		    int *n,
		    float *omega,float *plon,
		    float *plat,
		    float *rate)
{
  int i,j,m[3],iweight=0;
  float r1[3],r2[3],v1[3],v2[3],polar_base[3][3],
    theta,phi,a,b,c,d,e,f;
  const float RADIUS=6371000.0;

  *plat= *plon= *rate=0.0;
  for(i=0;i<3;i++){
    omega[i]=0.0;
    m[i]=0;
  }
  // need at least two points
  if(*n<2)
    return;
  // obtain r1 (cartesian)
  p2c(RADIUS,(theta=LAT2THETA(lat[0])),(phi=LON2PHI(lon[0])),r1);
  // obtain base vector
  calc_polar_base_at_r(polar_base,theta,phi);
  // obtain v1
  polar_vec2cart_vec(velt[0],velp[0],v1,polar_base);
  
  for(i=1;i< *n - 1;i++){
    // do the same for r2, v2
    p2c(RADIUS,(theta=LAT2THETA(lat[i])),
	(phi=LON2PHI(lon[i])),r2);
    calc_polar_base_at_r(polar_base,theta,phi);
    polar_vec2cart_vec(velt[i],velp[i],v2,polar_base);

    a=(r1[Z]*r2[X] - r1[X]*r2[Z]);
    b=(r1[Y]*r2[X] - r1[X]*r2[Y]);
    if(a!=0){
      omega[X] += -((  r2[X]*v1[Y] -  r1[X]*v2[Y])/a);
      m[X]++;
    }
    if(b!=0){
      omega[X] += -((-(r2[X]*v1[Z]) + r1[X]*v2[Z])/b);
      m[X]++;
    }

    c=(r1[Z]*r2[Y] - r1[Y]*r2[Z]);
    d=(r1[Y]*r2[X] - r1[X]*r2[Y]);
    if(c!=0){
      omega[Y] += -((-(r2[Y]*v1[X]) + r1[Y]*v2[X])/c);
      m[Y]++;
    }
    if(d!=0){
      omega[Y] += -((-(r2[Y]*v1[Z]) + r1[Y]*v2[Z])/d);
      m[Y]++;
    }

    e=(r1[Z]*r2[Y] - r1[Y]*r2[Z]);
    f=(r1[Z]*r2[X] - r1[X]*r2[Z]);
    if(e!=0){
      omega[Z] += -((-(r2[Z]*v1[X]) + r1[Z]*v2[X])/e);
      m[Z]++;
    }
    if(f!=0){
      omega[Z] += -((  r2[Z]*v1[Y]  - r1[Z]*v2[Y])/f);
      m[Z]++;
    }
    for(j=0;j<3;j++){
      r1[j]=r2[j];
      v1[j]=v2[j];
    }
  }
  for(i=0;i<3;i++){
    if(m[i])
      omega[i]/=(float)m[i];
    fprintf(stderr,"%6i/%6i\n",m[i],*n);
  }  
  xyz2rtp(omega,v1);
  *plon= PHI2LON(v1[PHI]);
  *plat= THETA2LAT(v1[THETA]);
  *rate= v1[R];


}
void xyz2rtp(float *x, float *r)
{
  double tmp1,tmp2;
  tmp1=SQUARE((double)x[X]) + SQUARE((double)x[Y]);
  tmp2=tmp1 + SQUARE((double)x[Z]);
  
  r[R]=(float)(sqrtl(tmp2));
  r[THETA]=(float)(atan2(sqrt(tmp1),(double)x[Z]));
  r[PHI]=(float)(atan2((double)x[Y],(double)x[X]));
  if(r[PHI]<0.0)
    r[PHI] += TWOPI;
  
}

void p2c(float r, 
	 float theta,
	 float phi,
	 float *xc)
{
  float tmpdbl;
  tmpdbl=sin(theta)*r;
  xc[X]=tmpdbl * cos(phi);
  xc[Y]=tmpdbl * sin(phi);
  xc[Z]=cos(theta)*r;
}

void polar_vec2cart_vec(float vtheta,float vphi,
			float *cart_vec,
			float polar_base[3][3])
{
  int i;
  for(i=0;i<3;i++){
    cart_vec[i]  = polar_base[THETA][i]* vtheta;
    cart_vec[i] += polar_base[PHI][i]  * vphi;
  }
}

void calc_polar_base_at_r(float polar_base[3][3],
			  float theta,
			  float phi)
{
  double ct,cp,st,sp;
  
  ct=cos(theta);
  cp=cos(phi);
  st=sin(theta);
  sp=sin(phi);
  
  polar_base[R][X]= st * cp;
  polar_base[R][Y]= st * sp;
  polar_base[R][Z]= ct;
  
  polar_base[THETA][X]= ct * cp;
  polar_base[THETA][Y]= ct * sp;
  polar_base[THETA][Z]= -st;
  
  polar_base[PHI][X]= -sp;
  polar_base[PHI][Y]= cp;
  polar_base[PHI][Z]= 0.0;
  


}

