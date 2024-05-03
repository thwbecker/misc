#include <math.h>
#include <stdio.h>

/* 

grid on an even area parameterization of a half-sphere

*/
#define PREC double
#define PIF 57.295779513082320876798154814105

struct r{
  int nc;
  PREC *lon;
  PREC dlon,lat;
  int *hit;
};
int main(int argc, char **argv)
{

  struct r *rows;
  PREC dlon0,readlon,readlat,x0,y0,x1,y1;
  int nr,nr0,i,j,nc0,ilon,ilat,nboxtotal,filled;
  
  dlon0 = 1;			/* start spacing on lon in degrees */
  
  nr = 90 / dlon0 ;		/* number of rows */
  nc0 = 360/dlon0;		/* number of columns at equator */
  
  /* get row array */
  rows = (struct r *)malloc(sizeof(struct r)*nr);
  /* intialize rows */
  nboxtotal = 0;
  for(i=0;i < nr;i++){
    rows[i].lat = (PREC)i * dlon0;
    rows[i].dlon = dlon0/cos(rows[i].lat/PIF);
    rows[i].nc = (PREC)360 / rows[i].dlon ;
    rows[i].lon = (PREC *)malloc(sizeof(PREC)*rows[i].nc);
    rows[i].hit = (int *)calloc(rows[i].nc,sizeof(int));
    nboxtotal += rows[i].nc;
    //fprintf(stderr,"%i %i %i %g %g\n",i,rows[i].nc,nc0,rows[i].lat,rows[i].dlon);
    for(j=0;j < rows[i].nc;j++){
      rows[i].lon[j] = rows[i].dlon * j;
      //fprintf(stdout,"%i(%i) %i(%i) lon: %g lat: %g\n",i,nr,j,rows[i].nc,rows[i].lon[j],rows[i].lat);
    }
  }
  /* initialization done */
  while(fscanf(stdin,"%lf %lf",&readlon,&readlat)==2){ /* read lon lat from stdin until done */
    if(readlon>360)
      readlon -= 360;
    if(readlon<0)
      readlon += 360;
    /* upper hemisphere */
    if(readlat < 0){
      fprintf(stderr,"%s: error: can only deal with positive latitudes!\n",argv[0]);
      exit(-1);
    }
    
    /* find row */
    ilat  =  (int)(readlat / dlon0);
    if(ilat >= nr){fprintf(stderr,"%s: error: out of bounds\n",argv[0]);exit(-1);}
    /* find column */
    ilon=0;
    while((ilon < rows[ilat].nc)&&(readlon > rows[ilat].lon[ilon]))
      ilon++;
    if(ilon== rows[ilat].nc)
      ilon=rows[ilat].nc;
    ilon--;
    rows[ilat].hit[ilon]++;	/* increment counter */

  } /* end of input loop */
  

  /* count how many boxes have entries */
  filled = 0;
  for(i=0;i<nr;i++)
    for(j=0;j<rows[i].nc;j++)
      if(rows[i].hit[j]){
	filled++;
	x0 = rows[i].lon[j];y0=rows[i].lat;
	if(j<rows[i].nc-1)
	  x1 = rows[i].lon[j+1];
	else
	  x1 = 360;
	if(i<nr-1)
	  y1 = rows[i+1].lat;
	else
	  y1 = 90;
	//fprintf(stdout,"%g %g\n%g %g\n%g %g\n%g %g\n>\n",
	//x0,y0,x1,y0,x1,y1,x0,y1);
      }
  //fprintf(stdout,"%s: %i out of %i boxes filled, ratio: %g\n",
  //  argv[0],filled,nboxtotal,(PREC)filled/(PREC)nboxtotal);
  printf("%g\n",(PREC)filled/(PREC)nboxtotal);
}
