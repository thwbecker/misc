#include <stdio.h>
#include <limits.h>
#include <math.h>
#include "determine_netr.h"

int main(int argc,char **argv)
{
  int n,zero=0;
  float *lon,*lat,plon,plat,rotvec[3],*velt,*velp,rate;
 
  switch(argc){
  case 1:{
    break;
  }
  default:{
    fprintf(stderr,"%s\nreads in lon[deg] lat[deg] v_x v_y from stdin\n",
	    argv[0]);
    fprintf(stderr,"\toutput is rotation pole in w_x w_y w_z\n");
    exit(-1);
    break;
  }}
  read_velocities(&lon,&lat,&velt,&velp,&n,0,argv[0]);
  fprintf(stderr,"%s: read in %i lon[deg] lat[deg] v_x v_y pairs\n",
	  argv[0],n);

#ifdef LINUX_SUBROUTINE_CONVENTION
#ifdef GCC_USCR
  determine_netr__(lon,lat,velt,velp,&n,rotvec,&plon,&plat,
		 &rate,&zero);
#else
  determine_netr_(lon,lat,velt,velp,&n,rotvec,&plon,&plat,
		  &rate,&zero);
#endif

#else  
  determine_netr(lon,lat,velt,velp,&n,rotvec,&plon,&plat,
		 &rate,&zero);
#endif
  fprintf(stderr,"%s: cartesian rotation pole:\n",argv[0]);
  fprintf(stdout,"%g %g %g\n",
	  rotvec[0],
	  rotvec[1],
	  rotvec[2]);

  return 0;
}
