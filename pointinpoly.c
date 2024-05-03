/* 

point in spherical polygon test modified from geometry gems 5



THIS DOES NOT WORK, WHO KNOWS WHY

*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "polysort.h"


int main ( int argc, char **argv ) 
{
 FILE     *f;
 int      nv,n,nin;
 struct GeoPoint *verts, p;
 int inside;
 Rdouble  Area,lon,lat;
 if(argc < 2){
   fprintf(stderr,"%s polyfile.dat\nreads lon lat in deg from stdin and determines if in polygon\n",
	   argv[0]);
   exit(-1);
 }

 /* read polygon file */
 if ( (f = fopen ( argv[1], "r" )) == NULL ){
   fprintf ( stdout, "%s: can not open the polyhedron file %s\n",
	     argv[0],argv[1]);
   exit ( -1 );
 }
 nv=0;
 verts = (struct GeoPoint *)malloc(sizeof(struct GeoPoint));
 while ( fscanf ( f, "%lf %lf", &lon,&lat) == 2 ){
   ll2xyz(lon,lat,(verts+nv));
   nv++;
   verts = (struct GeoPoint *)realloc(verts,sizeof(struct GeoPoint)*(nv+1));
 }
 fclose(f);
 fprintf(stderr,"%s: read %i nodes for polygon\n",argv[0],nv);
 /* done reading poly, start checking nodes */

 n=0;nin=0;
 while (fscanf( stdin, "%lf %lf", &lon,&lat ) == 2){
   ll2xyz(lon,lat,&p);
   inside = pointinpoly(nv, verts, &p);
   fprintf ( stdout, "%g %g %i\n", lon,lat,inside);
   n++;
   nin+=inside;
 }
 fprintf(stderr,"%s: %i out of %i nodes tested OK\n",
	 argv[0],nin,n);
 return 0;
}
