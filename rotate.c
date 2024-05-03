#include <math.h>
#include "rotate.h"
#include "coordconvention.h"
#include "trig.h"
/*

  given the azimuthal (CW from north) angle alpha (degrees), generate
  the rotation matrix that transforms between the coordinate systems

*/
void get_z_rot(COMP_PRECISION r[3][3],COMP_PRECISION alpha)
{
  COMP_PRECISION ralpha,sin_alpha,cos_alpha;

  ralpha=DEG2RAD(alpha);
  sin_alpha = sin(ralpha);
  cos_alpha = cos(ralpha);

  r[X][X] = cos_alpha;
  r[X][Y] = sin_alpha;
  r[X][Z] = 0.0;

  r[Y][X] = -sin_alpha;
  r[Y][Y] = cos_alpha;
  r[Y][Z] = 0.0;
  
  r[Z][X] = 0.0;
  r[Z][Y] = 0.0;
  r[Z][Z] = 1.0;
}
/*
  
  obtain a general rotation matrix with angles 
  alpha, beta, and gamma (given in degrees) as 
  defined in Dahlen and Tromp, p. 921
  
*/
void get_gen_rot(COMP_PRECISION r[3][3],COMP_PRECISION alpha,
		 COMP_PRECISION beta, COMP_PRECISION gamma)
{
  COMP_PRECISION ralpha,sin_alpha,cos_alpha;
  COMP_PRECISION rbeta,sin_beta,cos_beta;
  COMP_PRECISION rgamma,sin_gamma,cos_gamma;

  ralpha=DEG2RAD(alpha);
  sin_alpha = sin(ralpha);
  cos_alpha = cos(ralpha);

  rbeta=DEG2RAD(beta);
  sin_beta = sin(rbeta);
  cos_beta = cos(rbeta);

  rgamma=DEG2RAD(gamma);
  sin_gamma = sin(rgamma);
  cos_gamma = cos(rgamma);
  
  r[X][X] = cos_alpha*cos_beta*cos_gamma - sin_alpha*sin_gamma; 
  r[X][Y] = sin_alpha*cos_beta*cos_gamma + cos_alpha*sin_gamma;
  r[X][Z] = -sin_beta*cos_gamma;
  r[Y][X] = -cos_alpha*cos_beta*sin_gamma -sin_alpha*cos_gamma;
  r[Y][Y] = -sin_alpha*cos_beta*sin_gamma +cos_alpha*cos_gamma;
  r[Y][Z] = sin_beta*sin_gamma;
  r[Z][X] = cos_alpha*sin_beta;
  r[Z][Y] = sin_alpha*sin_beta;
  r[Z][Z] = cos_beta;
}

/*

  rotates a 3D vector (such as a displacement vector) into a new
  coordinate system whose x-axes is rotated around the 
  z-axis by alpha angle in clockwise direction from east

  WARNING: this routine doesn't take azimuth but direction!

  alternatively: rotates a vector in the old x-y-z coordinate system
  into a new one whose x-axis is rotated by an angle alpha from
  west-east (old x-axis) in a counterclockwise direction

*/
void rotate_vec_2d_angle(COMP_PRECISION *xin, COMP_PRECISION *xout,
			 COMP_PRECISION cos_alpha, COMP_PRECISION sin_alpha)
{
  xout[0]= ( cos_alpha * xin[0]) + (sin_alpha * xin[1]);
  xout[1]= (-sin_alpha * xin[0]) + (cos_alpha * xin[1]);
}
/* 
   rotate a 3-D vector
*/
void rotate_vec_3d(COMP_PRECISION xin[3],COMP_PRECISION xout[3],
		   COMP_PRECISION r[3][3])
{
  int i,j,k;
  // calculate r. xin
  for(i=0;i<3;i++){
    xout[i] = 0.0;
    for(j=0;j<3;j++){
      xout[i] += xin[j] * r[i][j];
    }
  }
}

/* 
   rotate a 3x3 tensor using a general rotation matrix r 
   xout = r . xin . r^T
*/
void rotate_mat(COMP_PRECISION xin[3][3],COMP_PRECISION xout[3][3],
		COMP_PRECISION r[3][3])
{
  COMP_PRECISION tmp[3][3];
  int i,j,k;
  // calculate xin . r^T
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      tmp[i][j] = 0.0;
      for(k=0;k<3;k++)
	tmp[i][j] += xin[i][k] * r[j][k];
    }
  // calculate r . tmp
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      xout[i][j] = 0.0;
      for(k=0;k<3;k++)
	xout[i][j] += r[i][k] * tmp[k][j];
    }
}
/* 
   
converts a a[3][3] matrix in polar system to a matrix b[3][3] in a cartesian
system, xp are the coordinates in the spherical system

*/
void polar_to_cart_mat_at_r3x3(COMP_PRECISION a[3][3],COMP_PRECISION b[3][3],
			       COMP_PRECISION *xp)
{
  COMP_PRECISION rot[3][3];
  /* obtain the rotation matrix at the location in spherical sys */
  calc_polar_basis_at_r(rot,xp); 
  rotate_mat(a,b,rot);	/* rotate tensor */
}
/* 

determine the polar basis vectors in cartesian space at position
r(r,theta,phi). this matrix has the r, theta, phi vectors in cart. components 
as rows R = ( e_r e_theta e_phi )

*/
void calc_polar_basis_at_r(COMP_PRECISION polar_basis[3][3],
			   COMP_PRECISION *r)
{
  COMP_PRECISION ct,cp,st,sp;

  my_sincos(r[PHI],&sp,&cp);
  my_sincos(r[THETA],&st,&ct);

  polar_basis[X][R]= st * cp;
  polar_basis[X][THETA]= ct * cp;
  polar_basis[X][PHI]= -sp;
  polar_basis[Y][R]= st * sp;
  polar_basis[Y][THETA]= ct * sp;
  polar_basis[Y][PHI]= cp;
  polar_basis[Z][R]= ct;
  polar_basis[Z][THETA]= -st;
  polar_basis[Z][PHI]= 0.0;
}
/* 
   
compute sines and cosines of x, given in radians

*/
void my_sincos(COMP_PRECISION x, COMP_PRECISION *sinv, COMP_PRECISION *cosv)
{
  *sinv = sin(x);
  *cosv = cos(x);
}


