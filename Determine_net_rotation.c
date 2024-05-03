//
//     calculates net rotation of input velocities
//     taken from code by Peter Bird, for copyright see end of routine
//
//     input is lon lat v_theta v_phi as n-dimensional vectors
//
#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>


void determine_netr_tp(float, float, float, float, float, int, double [3][3], double *);
void sub_netr(float, float, float, float *, float *, double *);
void hc_ludcmp_3x3(double [3][3], int *);
void hc_lubksb_3x3(double [3][3], int *, double *);


void main()
{
  int n,i;
  float *vphi,*vtheta,*r,*theta,*phi,vdummy;
  double coef[3][3],omega[3];
  

  vphi = (float *)malloc(sizeof(float));
  vtheta = (float *)malloc(sizeof(float));
  r = (float *)malloc(sizeof(float));
  theta = (float *)malloc(sizeof(float));
  phi = (float *)malloc(sizeof(float));

  n = 0;
  while(fscanf(stdin,"%f %f %f %f %f",(r+n),(theta+n),(phi+n),(vtheta+n),(vphi+n))==5){
    if(n == 0)
      determine_netr_tp(vdummy,vdummy,vdummy,vdummy,vdummy,0,coef,omega); /* init */

    determine_netr_tp(r[n],theta[n],phi[n],vtheta[n],vphi[n],1,coef,omega); /* add */
    n++;
    vphi = (float *)realloc(vphi,sizeof(float)*(n+1));
    vtheta = (float *)realloc(vtheta,sizeof(float)*(n+1));
    r = (float *)realloc(r,sizeof(float)*(n+1));
    theta = (float *)realloc(theta,sizeof(float)*(n+1));
    phi = (float *)realloc(phi,sizeof(float)*(n+1));
  }
  determine_netr_tp(vdummy,vdummy,vdummy,vdummy,vdummy,2,coef,omega); /* solve */

  fprintf(stderr,"read %i entries, net rotation: %g,%g,%g\n",n,omega[0],omega[1],omega[2]);
  /* remove */
  for(i=0;i<n;i++)
    sub_netr(r[i],theta[i],phi[i],(vtheta+i),(vphi+i),omega);
  /* compute new */
  determine_netr_tp(vdummy,vdummy,vdummy,vdummy,vdummy,0,coef,omega);
  for(i=0;i<n;i++)
    determine_netr_tp(r[i],theta[i],phi[i],vtheta[i],vphi[i],1,coef,omega); /* add */
  determine_netr_tp(vdummy,vdummy,vdummy,vdummy,vdummy,2,coef,omega); /* solve */
  fprintf(stderr,"after correction: net rotation: %g,%g,%g\n",omega[0],omega[1],omega[2]);


  
}

double determine_netr_tp(float r,float theta,float phi,
			 float velt,float velp,int mode,
			 double coef[3][3],double *omega)
{
  float coslat,coslon,sinlat,sinlon,rx,ry,rz,rate,rzu,a,b,c,d,e,f;
  int i,j,ind[3];
  
  switch(mode){
  case 0:			/* initialize */
    for(i=0;i < 3;i++){
      for(j=0;j < 3;j++)
	coef[i][j] = 0.0;
      omega[i] = 0.0;
    }
    break;
  case 1:			/* add this velocity */
    if((fabs(theta) > 1e-6) &&(fabs(theta-M_PI)>1e-6)){
      coslat=sin(theta);
      coslon=cos(phi);
      sinlat=cos(theta);
      sinlon=sin(phi);
      
      rx=coslat*coslon*r;
      ry=coslat*sinlon*r;
      rz=sinlat*r;
      
      rzu=sinlat;

      a = -rz*rzu*sinlon-ry*coslat;
      b = -rz*coslon;
      c =  rz*rzu*coslon+rx*coslat;
      d = -rz*sinlon;
      e = -ry*rzu*coslon+rx*rzu*sinlon;
      f =  ry*sinlon+rx*coslon;
      
      coef[0][0] += a*a+b*b;
      coef[0][1] += a*c+b*d;
      coef[0][2] += a*e+b*f;
      coef[1][1] += c*c+d*d;
      coef[1][2] += c*e+d*f;
      coef[2][2] += e*e+f*f;
      
      omega[0] += a*velt+b*velp;
      omega[1] += c*velt+d*velp;
      omega[2] += e*velt+f*velp;
    }
    break;
  case 2:			/* solve */
    coef[1][0]=coef[0][1];
    coef[2][0]=coef[0][2];
    coef[2][1]=coef[1][2];
    /* solve and overwrite omega with solution*/
    hc_ludcmp_3x3(coef,ind);
    hc_lubksb_3x3(coef,ind,omega);
    break;
  default:
    fprintf(stderr,"determine_netr_tp: mode %i undefined\n",mode);
    exit(-1);
    break;
  }
}

//     subtract a net rotation component from a velocity 
//     field
      
void sub_netr(float r,float theta,float phi,float *velt,float *velp, double *omega)
{

  float coslat,coslon,sinlon,sinlat,rx,ry,rz,rzu;
  float vx,vy,vz,tx,ty,tz,vtheta,pc,px,py,vphi;

  coslat=sin(theta);
  coslon=cos(phi);
  sinlat=cos(theta);
  sinlon=sin(phi);

  rx=coslat*coslon*r;
  ry=coslat*sinlon*r;
  rz=sinlat*r;

  rzu=sinlat;
  
  vx = omega[1]*rz - omega[2]*ry;
  vy = omega[2]*rx - omega[0]*rz;
  vz = omega[0]*ry - omega[1]*rx;
  
  tx=rzu*coslon;
  ty=rzu*sinlon;
  tz=-coslat;

  vtheta=vx*tx+vy*ty+vz*tz;

  px=-sinlon;
  py=coslon;
  vphi = vx*px+vy*py;

  *velt -= vtheta;
  *velp -= vphi;
}


  
//
//      PROGRAM -OrbScore-: COMPARES OUTPUT FROM -SHELLS-
//                          WITH DATA FROM GEODETI// NETWORKS,
//                          STRESS DIRECTIONS, FAULT SLIP RATES,
//                          SEAFLOOR SPREADING RATES, AND SEISMICITY,
//                          AND REPORTS SUMMARY SCALAR SCORES.
//
//=========== PART OF THE "SHELLS" PACKAGE OF PROGRAMS===========
//
//   GIVEN A FINITE ELEMENT GRID FILE, IN THE FORMAT PRODUCED BY
//  -OrbWeave- AND RENUMBERED BY -OrbNumbr-, WITH NODAL DATA
//   ADDED BY -OrbData-, AND NODE-VELOCITY OUTPUT FROM -SHELLS-,
//   COMPUTES A VARIETY OF SCORES OF THE RESULTS.
//
//   NOTE: Does not contain VISCOS or DIAMND, hence independent
//         of changes made in May 1998, and equally compatible
//         with Old_SHELLS or with improved SHELLS.
//
//                             by
//                        Peter Bird
//          Department of Earth and Spcae Sciences,
//    University of California, Los Angeles, California 90095-1567
//   (C) Copyright 1994, 1998, 1999, 2000
//                by Peter Bird and the Regents of
//                 the University of California.
//              (For version data see FORMAT 1 below)
//
//   THIS PROGRAM WAS DEVELOPED WITH SUPPORT FROM THE UNIVERSITY OF
//     CALIFORNIA, THE UNITED STATES GEOLOGI// SURVEY, THE NATIONAL
//     SCIENCE FOUNDATION, AND THE NATIONAL AERONAUTICS AND SPACE
//     ADMINISTRATION.
//   IT IS FREEWARE, AND MAY BE COPIED AND USED WITHOUT CHARGE.
//   IT MAY NOT BE MODIFIED IN A WAY WHICH HIDES ITS ORIGIN
//     OR REMOVES THIS MESSAGE OR THE COPYRIGHT MESSAGE.
//   IT MAY NOT BE RESOLD FOR MORE THAN THE COST OF REPRODUCTION
//      AND MAILING.
//     



/* 

matrix solvers from numerical recipes

 */
#define NR_TINY 1.0e-20;

void hc_ludcmp_3x3(double a[3][3],int *indx)
{
  int i,imax=0,j,k;
  double big,dum,sum,temp;
  double vv[3];
  
  for (i=0;i < 3;i++) {
    big=0.0;
    for (j=0;j < 3;j++)
      if ((temp = fabs(a[i][j])) > big) 
	big=temp;
    if (fabs(big) < 5e-15) {
      fprintf(stderr,"hc_ludcmp_3x3: singular matrix in routine, big: %g\n",
	      big);
      //hc_print_3x3(a,stderr);
      for(j=0;j<3;j++)
	fprintf(stderr,"%g %g %g\n",a[j][0],a[j][1],a[j][2]);
      exit(-1);
    }
    vv[i]=1.0/big;
  }
  for (j=0;j < 3;j++) {
    for (i=0;i < j;i++) {
      sum = a[i][j];
      for (k=0;k < i;k++) 
	sum -= a[i][k] * a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i < 3;i++) {
      sum=a[i][j];
      for (k=0;k < j;k++)
	sum -= a[i][k] * a[k][j];
      a[i][j]=sum;
      if ( (dum = vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k < 3;k++) {
	dum = a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (fabs(a[j][j]) < 5e-15) 
      a[j][j] = NR_TINY;
    if (j != 2) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i < 3;i++) 
	a[i][j] *= dum;
    }
  }
}
#undef NR_TINY
void hc_lubksb_3x3(double a[3][3], int *indx, double *b)
{
  int i,ii=0,ip,j;
  double sum;
  for (i=0;i < 3;i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii-1;j <= i-1;j++) 
	sum -= a[i][j]*b[j];
    else if (fabs(sum) > 5e-15) 
      ii = i+1;
    b[i]=sum;
  }
  for (i=2;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j < 3;j++) 
      sum -= a[i][j]*b[j];
    b[i] = sum/a[i][i];
  }
}

