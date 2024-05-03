#include "catalog.h"
/* 


perform a smoothed seismicity summation

mode 0: scaled with moment 
mode 1: scaled with sqrt(moment)
mode 2: scaled with unity

*/

int main(int argc, char **argv)
{
  double dx = 0.2, dy = 0.2,min_mag=2.49,max_mag=8.01;
  double lonmin,lonmax,latmin,latmax;
  double mindepth,maxdepth;
  struct cat *catalog;
  char outfile[300];
  int mode = 0;

  catalog=(struct cat *)malloc(sizeof(struct cat)); 

  if(!catalog)MEMERROR(argv[0]);

  /* 
     range of summation 
  */
  lonmin = 232;lonmax=250;
  latmin=30;latmax = 50;
  mindepth = 0;maxdepth = 15;	/* in km */
  
  if(argc < 2){
    fprintf(stderr,"%s catalog.xyzmt [mode, %i] [dx/dy, %g] [min_mag, %g] [max_mag, %g] [min_lon, %g] [max_lon, %g] [min_lat, %g] [max_lat, %g] [max_depth, %g]\n",
	    argv[0],
	    mode,dx,min_mag,max_mag,
	    lonmin, lonmax, 
	    latmin, latmax,
	    maxdepth);
    exit(-1);
  }
  if(argc>2){
    sscanf(argv[2],"%i",&mode);
  }
  if(argc>3){
    sscanf(argv[3],"%lf",&dx);
    dy = dx;
  }
  if(argc>4)
    sscanf(argv[4],"%lf",&min_mag);
  if(argc>5)
    sscanf(argv[5],"%lf",&max_mag);
  if(argc>6)sscanf(argv[6],"%lf",&lonmin);
  if(argc>7)sscanf(argv[7],"%lf",&lonmax);
  if(argc>8)sscanf(argv[8],"%lf",&latmin);
  if(argc>9)sscanf(argv[9],"%lf",&latmax);
  if(argc>10)sscanf(argv[10],"%lf",&maxdepth);

  sprintf(outfile,"sseis.%g.%i",dx,mode);

  /* 
     read events
   */
  read_catalog(argv[1],catalog,ENG);

  /* 
     
     setup bins
     
  */
  /* parameters */
  /*  */
  setup_kostrov(catalog,lonmin,lonmax,latmin,latmax,mindepth,maxdepth,
		dx,dy,min_mag,max_mag);  
  /* 
     sum 
  */
  sum_smoothed_seismicity(catalog,mode);
  /* 
     print non-zero summed entries
  */
  print_summed_moment(catalog, outfile);
  return 0;
}

