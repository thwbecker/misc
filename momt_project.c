#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/* 

   reads a CMT moment tensor in 

   
   X Y depth mrr mtt mff mrt mrf mtf exp newX newY time_second 

   format and rotates it to pscoupe style side view along a profile with azimuth azi CW from NW

   


*/

#define X 0
#define Y 1
#define Z 2

#define PIOVERONEEIGHTY 0.0174532925199433
#define DEG2RAD(x) ( (x)*PIOVERONEEIGHTY )

#define COMP_PRECISION double

void rotate_mat(COMP_PRECISION [3][3],COMP_PRECISION [3][3],
		COMP_PRECISION [3][3]);
void get_z_rot(COMP_PRECISION [3][3],COMP_PRECISION );
void get_gen_rot(COMP_PRECISION [3][3],COMP_PRECISION,
		 COMP_PRECISION,COMP_PRECISION );


int main(int argc, char **argv)
{
  COMP_PRECISION azi,r[3][3],xin[3][3],xout[3][3];

  COMP_PRECISION x, y, depth, newx, newy, time_second;
  int exp;
  if((argc<2)||(argc>4)){
    fprintf(stderr,"%s azi[deg]\n\n",argv[0]);
    fprintf(stderr,"reads a CMT moment tensor in\n\n");
    fprintf(stderr,"X Y depth mrr mtt mff mrt mrf mtf exp newX newY time_second\n\n");
    fprintf(stderr,"format and rotates it to pscoupe style side view along a profile with azimuth azi CW from NW\n");
    exit(-1);
  }
  sscanf(argv[1],"%lf",&azi);
  
  /* make the rotation matrix to view from side with depth */
  get_gen_rot(r,azi,-90,90);

  //
  // X Y depth
  // mrr mtt mff
  // mrt mrf mtf
  // exp newX newY time_second 
  while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %i %lf %lf %lf\n",
	       &x,&y,&depth,
	       &xin[Z][Z],&xin[Y][Y],&xin[X][X],
	       &xin[Y][Z],&xin[X][Z],&xin[X][Y],
	       &exp,&newx,&newy,&time_second)==13){
    // flip the sign of the components with theta, since t = -y
    xin[Y][Z] *= -1.0;
    xin[X][Y] *= -1.0;
    /* fill in */
    xin[Y][X]=xin[X][Y];xin[Z][X]=xin[X][Z];xin[Z][Y]=xin[Y][Z];

    rotate_mat(xin,xout,r);
    /* output */
    fprintf(stdout,"%11g %11g %11g %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %3i %11g %13g %9i\n",
	    x,y,depth,
	    xout[Z][Z],xout[Y][Y],xout[X][X],
	    -xout[Z][Y],xout[Z][X],-xout[Y][X],
	    exp,newx,newy,(int)time_second);

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
