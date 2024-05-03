#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//
// converts data in lon lat v_r v_theta v_phi 
// format from polar vector components to cartesian components
// input is lon lat v_r v_theta v_phi
// output is lon lat v_x v_y v_z
// $Id: lonlatpv2cv.c,v 1.4 2001/10/17 22:35:19 becker Exp $
//
#define PI 3.1415926535897932384626434

int main(int argc, char **argv)
{
  double f,theta,phi,vr,vtheta,vphi,ct,st,cp,sp,
    polar_base_r[3],polar_base_theta[3],polar_base_phi[3],
    cart_vec[3],lon,lat;
  int i;

  f=PI/180.0;
  
  if(argc !=1 ){
    fprintf(stderr,"reads lon lat v_r v_theta v_phi and writes lon lat v_x v_y v_z\n");
    exit(-1);
  }
  

  while(fscanf(stdin,"%lf %lf %lf %lf %lf",&lon,&lat,&vr,&vtheta,&vphi)==5){

    theta=(90.0-lat)*f;
    phi=lon*f;

    // base vecs
    ct=cos(theta);cp=cos(phi);
    st=sin(theta);sp=sin(phi);
    //
    polar_base_r[0]= st * cp;
    polar_base_r[1]= st * sp;
    polar_base_r[2]= ct;
    //
    polar_base_theta[0]= ct * cp;
    polar_base_theta[1]= ct * sp;
    polar_base_theta[2]= -st;
    //
    polar_base_phi[0]= -sp;
    polar_base_phi[1]=  cp;
    polar_base_phi[2]= 0.0;
    // convert vector
    for(i=0;i<3;i++){
      cart_vec[i]  = polar_base_r[i]    * vr;
      cart_vec[i] += polar_base_theta[i]* vtheta;
      cart_vec[i] += polar_base_phi[i]  * vphi;
    }
    fprintf(stdout,"%20.16e %20.16e %20.16e %20.16e %20.16e\n",
	    lon,lat,cart_vec[0],cart_vec[1],cart_vec[2]);
  }
  return 0;
}








