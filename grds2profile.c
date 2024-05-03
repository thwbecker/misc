#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ggrd_grdtrack_util.h"
/* 

   read a set of scalar netcdf grids and provide profiles using the
   ggrd package

*/
#define DZ 1
int main(int argc, char **argv)
{
  struct ggrd_gt grids[1];
  ggrd_boolean verbose = FALSE;
  ggrd_boolean change_z_sign = TRUE;
  ggrd_boolean use_nearneighbor = FALSE;
  double lon,lat,z,value;
  if(argc < 3){
    fprintf(stderr,"%s - reads sets of Netcdf grd files and computes a profile\n",argv[0]);
    fprintf(stderr,"%s directory/dv directory/depth.dat\n",argv[0]);
    exit(-1);
  }
  if(ggrd_grdtrack_init_general(TRUE,argv[1],argv[2],"-fg",grids,
				 verbose,change_z_sign,use_nearneighbor)){
    fprintf(stderr,"%s: grd init error, grd: %s dfile: %s \n",argv[0],argv[1],argv[2]);
    exit(-2);
  }
  while(fscanf(stdin,"%lf %lf",&lon,&lat)==2){
    for(z=2871;z>=DZ;z-=DZ){
      ggrd_grdtrack_interpolate_lonlatz(lon,lat,z,grids,&value,verbose);
      printf("%11g %11g\t",z,value);
    }
    printf("\n");
  }

  return 0;
}
