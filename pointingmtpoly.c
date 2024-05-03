/* 

point in GMT polygons, uses GMT4

*/
#include "polysort.h"


int main ( int argc, char **argv ) 
{
  int n,nin,i,foundone;
  int inside;
  BOOLEAN greenwich = FALSE;
  double lon,lat;
  struct GMT_TABLE *pol;
  if(argc < 2){
    fprintf(stderr,"%s polyfile.gmt\nreads lon lat in deg from stdin and determines if in GMT style polygons in polyfile.gmt\n",
	    argv[0]);
    exit(-1);
  }
  argc = GMT_begin (argc, argv);
  GMT_io.in_col_type[0] =  GMT_io.out_col_type[0] = GMT_IS_LON;
  GMT_io.in_col_type[1] =  GMT_io.out_col_type[1] = GMT_IS_LAT;


  GMT_import_table (argv[1], GMT_IS_FILE, &pol, 0.0, greenwich, TRUE, FALSE);
  fprintf(stderr,"%s: read %i polygons\n",argv[0],pol->n_segments);
  /* done reading poly, start checking nodes */

  n=0;nin=0;
  while (fscanf( stdin, "%lf %lf", &lon,&lat ) == 2){
    foundone = 0;
    for(i=0;i < pol->n_segments;i++){
      inside = gmt_point_in_poly(pol,i,lon,lat);
      if(inside){		/* point is in this polygon */
	foundone  = 1;
	fprintf ( stdout, "%g %g %i\n", lon,lat,i+1);
	nin ++;
      }
    } /* end polygon loop */
    if(!foundone)
      fprintf ( stdout, "%g %g %i\n", lon,lat,0);
    n++;
  }
  fprintf(stderr,"%s: %i out of %i nodes tested inside polygons\n",
	  argv[0],nin,n);

  GMT_end (argc, argv);
  return 0;
}
