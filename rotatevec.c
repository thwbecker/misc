#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rotate.h"

/* 

   rotate a vector 3-D A by angles alpha and beta and gamma
   alpha (CCW from East (x-y), in degrees) around z or r axis, 
   (see dahlen and tromp, page 921)

   input:

   A_x A_y A_z for CARTESIAN

   A_r A_t A_p  for SPHERICAL

   output:

   A_x' A_y' A_z'   for CARTESIAN
   A_r' A_t' A_p'   for SPHERICAL

   where t and p are theta and phi in the original system
   and t' and p' are theta and phi in the rotated system
   for x and y the same



   REMEMBER THAT R^1(R(alpha,beta,gamma) is R(-gamma,-beta,-alpha)

   $Id: rotatevec.c,v 1.4 2005/02/23 03:04:21 becker Exp becker $

*/
#include "coordconvention.h"
#include "trig.h"
#include "rotate.h"

int main(int argc, char **argv)
{
  COMP_PRECISION alpha,beta,gamma,r[3][3],xin[3],xout[3];
  if((argc<2)||(argc>4) ){
    fprintf(stderr,"%s alpha[deg] [beta, deg] [gamma, deg]\n",
	    argv[0]);
#ifdef SPHERICAL
    fprintf(stderr,"rotate a 3-D vector A by angle alpha\n(CCW from East (x-y)) around r axis,\n");
    fprintf(stderr,"\n\tinput\n\tA_r A_t A_p \n\n\toutput:\n\n\t");
    fprintf(stderr,"A_r' A_t' A_p' \n");
    fprintf(stderr,"\nwhere t and p are theta and phi in the original system\n");
    fprintf(stderr,"and t' and p' are theta and phi in the rotated system\n\n");
#else// cartesian
    fprintf(stderr,"rotate a 3-D vector A by angle alpha\n(CCW from East (x-y)) around z axis,\n");
    fprintf(stderr,"and then by beta and gamma as in the Euler angle description of Dahlen & Tromp\n");
    fprintf(stderr,"\n\tinput\n\tA_x A_y A_z \n\n\toutput:\n\n\t");
    fprintf(stderr,"A_x' A_y' A_z'\n");
    fprintf(stderr,"\nwhere x, y, are z in the original system\n");
    fprintf(stderr,"and x', y', are z' in the rotated system\n\n");
#endif
    exit(-1);
  }
  /* 
     read angles in degrees
  */
  sscanf(argv[1],"%lf",&alpha);
  if(argc>2)
    sscanf(argv[2],"%lf",&beta);
  else
    beta=0.0;
  if(argc>3)
    sscanf(argv[3],"%lf",&gamma);
  else
    gamma=0.0;
  /* get the rotation matrix */
  get_gen_rot(r,alpha,beta,gamma);
  //#define VERBOSE
#ifdef VERBOSE
  fprintf(stderr,"%s: rotating system by alpha %g beta %g gamma %g degrees\n",
  	  argv[0],alpha,beta,gamma);
  fprintf(stderr,"%g %g %g\n %g %g %g\n %g %g %g\n",
	  r[X][X],r[X][Y],r[X][Z],r[Y][X],
	  r[Y][Y],r[Y][Z],r[Z][X],r[Z][Y],r[Z][Z]);
#endif
#ifdef SPHERICAL  
  while(fscanf(stdin,"%lf %lf %lf\n",&xin[Z],&xin[Y],&xin[X])==3){
    // flip the sign of the components with theta, since t = -y
    xin[Y] *= -1.0;
    rotate_vec_3d(xin,xout,r);
    fprintf(stdout,"%15.8e %15.8e %15.8e\n",
	    xout[Z],-xout[Y],xout[X]);
  }
#else
  while(fscanf(stdin,"%lf %lf %lf\n",&xin[X],&xin[Y],&xin[Z])==3){
    rotate_vec_3d(xin,xout,r);
    fprintf(stdout,"%15.8e %15.8e %15.8e\n",
	    xout[X],xout[Y],xout[Z]);
  }
#endif
  return 0;
}

