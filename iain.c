#include <math.h>
#include <stdio.h>


/* 

   rotate a symmetric 3x3 matrix A by angle alpha gamma
   (CCW from East (x-y), in degrees) around z or r axis, 
   see dahlen and tromp p. 921

   input:

   A_xx A_xy A_xz A_yy A_yz A_zz for CARTESIAN

   A_rr A_rt A_rp A_tt A_tp A_pp for SPHERICAL

   output:

   A_x'x' A_x'y' A_x'z A_y'y' A_y'z  A_zz   for CARTESIAN
   A_rr   A_rt'  A_rp' A_t't' A_t'p' A_p'p' for SPHERICAL

   where t and p are theta and phi in the original system
   and t' and p' are theta and phi in the rotated system
   for x and y the same


   REMEMBER THAT R^1(R(alpha,beta,gamma) is R(-gamma,-beta,-alpha)


   $Id: rotatemat.c,v 1.4 2005/02/23 03:04:21 becker Exp becker $

*/

#define X 0
#define Y 1
#define Z 2

#define PIOVERONEEIGHTY 0.0174532925199433
#define DEG2RAD(x) ( (x)*PIOVERONEEIGHTY )


void rotate_mat(COMP_PRECISION [3][3],COMP_PRECISION [3][3],
		COMP_PRECISION [3][3]);
void get_z_rot(COMP_PRECISION [3][3],COMP_PRECISION );
void get_gen_rot(COMP_PRECISION [3][3],COMP_PRECISION,
		 COMP_PRECISION,COMP_PRECISION );


int main(int argc, char **argv)
{
  COMP_PRECISION alpha,beta,gamma,r[3][3],xin[3][3],xout[3][3];
  if((argc<2)||(argc>4)){
    fprintf(stderr,"%s alpha[deg] [beta, deg] [gamma, deg]\n",argv[0]);
#ifdef SPHERICAL
    fprintf(stderr,"rotate a symmetric 3x3 matrix A by angle alpha\n(CCW from East (x-y)) around r axis,\n");
    fprintf(stderr,"\n\tinput\n\tA_rr A_rt A_rp A_tt A_tp A_pp\n\n\toutput:\n\n\t");
    fprintf(stderr,"A_rr A_rt' A_rp' A_t't' A_t'p' A_p'p'\n");
    fprintf(stderr,"\nwhere t and p are theta and phi in the original system\n");
    fprintf(stderr,"and t' and p' are theta and phi in the rotated system\n\n");
#else// cartesian
    fprintf(stderr,"rotate a symmetric 3x3 matrix A by angle alpha\n(CCW from East (x-y)) around z axis,\n");
    fprintf(stderr,"\n\tinput\n\tA_xx A_xy A_xz A_yy A_yz A_zz\n\n\toutput:\n\n\t");
    fprintf(stderr,"A_x'x' A_x'y' A_x'z A_y'y' A_y'z A_zz\n");
    fprintf(stderr,"\nwhere x and y are in the original system\n");
    fprintf(stderr,"and x' and y' in the rotated system\n\n");
#endif
    exit(-1);
  }
  sscanf(argv[1],"%lf",&alpha);
  if(argc>2)
    sscanf(argv[2],"%lf",&beta);
  else
    beta = 0.0;
  if(argc>3)
    sscanf(argv[3],"%lf",&gamma);
  else
    gamma = 0.0;
  /* 
     get the rotation matrix 
  */
  //fprintf(stderr,"%s: rotating system by alpha %g beta %g gamma %g degrees\n",
  //	  argv[0],alpha,beta,gamma);
  get_gen_rot(r,alpha,beta,gamma);
  //printf("%g %g %g %g %g %g %g %g %g\n",r[X][X],r[X][Y],r[X][Z],r[Y][X],r[Y][Y],r[Y][Z],r[Z][X],r[Z][Y],r[Z][Z]);
#ifdef SPHERICAL  
  while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf\n",
	       &xin[Z][Z],&xin[Y][Z],&xin[X][Z],
	       &xin[Y][Y],&xin[X][Y],&xin[X][X])==6){
    // flip the sign of the components with theta, since t = -y
    xin[Y][Z] *= -1.0;
    xin[X][Y] *= -1.0;
    xin[Y][X]=xin[X][Y];xin[Z][X]=xin[X][Z];xin[Z][Y]=xin[Y][Z];

    /* rotate the matrix */
    rotate_mat(xin,xout,r);
    fprintf(stdout,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
	    xout[Z][Z],-xout[Z][Y],xout[Z][X],
	       xout[Y][Y],-xout[Y][X],xout[X][X]);
#else
  while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf\n",
	       &xin[X][X],&xin[X][Y],&xin[X][Z],
	       &xin[Y][Y],&xin[Y][Z],&xin[Z][Z])==6){
    xin[Y][X]=xin[X][Y];xin[Z][X]=xin[X][Z];xin[Z][Y]=xin[Y][Z];
    rotate_mat(xin,xout,r);
    fprintf(stdout,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
	    xout[X][X],xout[X][Y],xout[X][Z],
	    xout[Y][Y],xout[Y][Z],xout[Z][Z]);
#endif
  }
  return 0;
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
