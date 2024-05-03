/* 


average splitting data in lon lat azi dt format 

according to GMT polygons

*/
#include "polysort.h"



int main ( int argc, char **argv ) 
{
  int      n,nin,i,nfound;
  int inside;
  BOOLEAN greenwich = FALSE;
  double lon,lat,azi,dt;
  struct GMT_TABLE *pol;
  struct poly *tri = NULL;

  if(argc < 2){
    fprintf(stderr,"%s poly.gmt\nreads lon lat azi[deg] dt  from stdin and\n",
	    argv[0]);
    fprintf(stderr,"\t\t\taverages splitting data in the GMT polygons  in poly.gmt file\n");
    exit(-1);
  }
  argc = GMT_begin (argc, argv);
  GMT_io.in_col_type[0] =  GMT_io.out_col_type[0] = GMT_IS_LON;
  GMT_io.in_col_type[1] =  GMT_io.out_col_type[1] = GMT_IS_LAT;


  GMT_import_table (argv[1], GMT_IS_FILE, &pol, 0.0, greenwich, TRUE, FALSE);
  fprintf(stderr,"%s: read %i GMT polygons\n",argv[0],pol->n_segments);
  /* init the tri structure, inlcuding centroids */
  gmt_poly_helper_init(pol,&tri);


  /* done reading poly, start checking nodes */
  
  n=nin=0;
  while (fscanf( stdin, "%lf %lf %lf %lf", 
		 &lon,&lat,&azi,&dt ) == 4){
    if((dt < 0)||(dt > 6)||(azi>360)){
      fprintf(stderr,"%s: node %i range error: azxi: %g (deg) dt: %g (s)\n",
	      argv[0],n,azi,dt);
      exit(-1);
    }
    /* 
       convert 
    */
    azi *= 2*RFAC;
    nfound=0;
    for(i=0;i < pol->n_segments;i++){		/* look for polygon */
      inside = gmt_point_in_poly(pol,i,lon,lat);
      if(inside){
	nfound++;
	nin++;			/* totoal nodes found */
	tri[i].azi = (Rdouble *)realloc(tri[i].azi,(tri[i].nin+1)*sizeof(Rdouble));
	tri[i].dt = (Rdouble *)realloc(tri[i].azi,(tri[i].nin+1)*sizeof(Rdouble));
	tri[i].azi[tri[i].nin] = azi;
	tri[i].dt[tri[i].nin] = dt;
	tri[i].nin++;		/* in this triangle */
	/* allow multiple assignments, e.g. when on border */
      }
    }
    if(!nfound){
      fprintf(stderr,"not associated: %g %g\n",lon,lat);
    }
    n++;
  }
  fprintf(stderr,"%s: processed  %i nodes, associated %i (allowing multiples)\n",
	  argv[0],n,nin);
  /* output */
  for(i=0;i < pol->n_segments;i++){
    //if(tri[i].nin){
      fprintf(stdout,"%11g %11g %i\n",
	      tri[i].lon_c,tri[i].lat_c,tri[i].nin);
      //    }
  }


  
  GMT_end (argc, argv);

  return 0;
}
