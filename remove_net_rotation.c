#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "determine_netr.h"

int main(int argc,char **argv)
{
  int n,use_weights=0,i;
  float *lon,*lat,plon,plat,rotvec[3],*velt,*velp,
    rate;
  switch(argc){
  case 1:{
    use_weights=0;
    break;}
  case 2:{
    sscanf(argv[1],"%i",&use_weights);
    break;}
  default:{
    fprintf(stderr,"%s [use_weights, 0]\n\t reads in lon[deg] lat[deg] v_theta v_phi from stdin\n",
	    argv[0]);
    fprintf(stderr,"\tdetermines the best fitting net rotation, and removes it\n");
    fprintf(stderr,"\toutput is lon[deg] lat[deg] v_theta v_phi to stdout\n");
    exit(-1);
    break;
  }}
  read_velocities(&lon,&lat,&velt,&velp, &n,use_weights,
		  argv[0]);
#ifdef LINUX_SUBROUTINE_CONVENTION
#ifdef GCC_USCR
  sub_netr__(lon,lat,velt,velp,&n,rotvec,&plon,&plat,&rate,&use_weights);
#else
  sub_netr_(lon,lat,velt,velp,&n,rotvec,&plon,&plat,&rate,&use_weights);
#endif
#else
  sub_netr(lon,lat,velt,velp,&n,rotvec,&plon,&plat,&rate,&use_weights);
#endif
  fprintf(stderr,"%s: best fitting net rotation: lon: %g lat: %g  magn: %g\n",
	  argv[0],plon,plat,rate);
  fprintf(stderr,"%s: wx %g wy: %g wz: %g\n",
	  argv[0],rotvec[0],rotvec[1],rotvec[2]);
  fprintf(stderr,"%s: output is reduced lon[deg] lat[deg] v_x v_y\n",argv[0]);
  for(i=0;i<n;i++)
    fprintf(stdout,"%g %g %g %g\n",lon[i],lat[i],velp[i],-velt[i]);

  return 0;
}
