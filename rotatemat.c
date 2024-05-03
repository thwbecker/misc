#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rotate.h"
#include <string.h>
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
#include "coordconvention.h"
#include "trig.h"
void rotate_mat(COMP_PRECISION [3][3],COMP_PRECISION [3][3],
		COMP_PRECISION [3][3]);

int main(int argc, char **argv)
{
  COMP_PRECISION alpha,beta,gamma,r[3][3],xin[3][3],xout[3][3];
  int alpha_from_file = 0;
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
    fprintf(stderr,"if alpha is set to \"file\" will expect one more column and read angles from file\n");
    exit(-1);
  }
  if(strcmp(argv[1],"file")==0)
    alpha_from_file = 1;
  else{
    sscanf(argv[1],"%lf",&alpha);
    alpha_from_file = 0;
  }
#ifdef SPHERICAL  
  if(alpha_from_file){
    fprintf(stderr,"alpha from file only implemented for cartesian\n");
    exit(-1);
  }
#endif
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
  if(!alpha_from_file)
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
  }
#else
  /* cartesian (what's with the spherical anyway?) */
  if(alpha_from_file){
    //fprintf(stderr,"%s: expecting alpha in seventh column\n",argv[0]);

    while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf %lf\n",
		 &xin[X][X],&xin[X][Y],&xin[X][Z],
		 &xin[Y][Y],&xin[Y][Z],&xin[Z][Z],&alpha)==7){
      get_gen_rot(r,alpha,beta,gamma);
      xin[Y][X]=xin[X][Y];xin[Z][X]=xin[X][Z];xin[Z][Y]=xin[Y][Z];
      rotate_mat(xin,xout,r);
      fprintf(stdout,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
	      xout[X][X],xout[X][Y],xout[X][Z],
	      xout[Y][Y],xout[Y][Z],xout[Z][Z]);
    }

  }else{
    while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf\n",
		 &xin[X][X],&xin[X][Y],&xin[X][Z],
		 &xin[Y][Y],&xin[Y][Z],&xin[Z][Z])==6){
      xin[Y][X]=xin[X][Y];xin[Z][X]=xin[X][Z];xin[Z][Y]=xin[Y][Z];
      rotate_mat(xin,xout,r);
      fprintf(stdout,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
	      xout[X][X],xout[X][Y],xout[X][Z],
	      xout[Y][Y],xout[Y][Z],xout[Z][Z]);
    }
  }
#endif
  return 0;
}

