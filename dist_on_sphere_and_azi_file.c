#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* 

compute distance on sphere and azimuth in degrees for two points
as read in from stdin, all in degrees

lon0 lat0 lon lat

output: 

dist azi lon lat

where dist and azimuth are in degrees, and the event location is merely repeated

*/
#define PIF 57.2957795130823
#define TWOPI 6.2831853071795864769252867665590
double distance(double ,double ,double ,double , double ,double);
#define stadis stadis_
extern void stadis(double *,double *,double *,double *,double *,double *);

int main(int argc, char **argv)
{
  double lon0,lat0,lon0d,lat0d,colat0d,lat,lon,slat[2],clat[2],dist,azi,lond,latd,colatd;
  static int mode = 1;		/* 0: my old version 1: new version using stadis */
  static short int verbose = 1;
  int azimode = 0,n;
  if(argc > 2){
    fprintf(stderr,"usage:\n\n%s [azi, 0]\n\ncompute distance and back-azimuth in degrees.\n",argv[0]);
    fprintf(stderr,"read lon0 lat0 lon lat from file\n");
    fprintf(stderr,"back-azimuth is the azimuth from the station to the source\n");
    fprintf(stderr,"if azi is set to unity, will compute azimuth from source to station instead\n");
    fprintf(stderr,"between lon0, lat0 (station) and a list of sources read from stdin, all in deg.\n");
    fprintf(stderr,"output: distance azimuth lon_source lat_source [deg]\n");
    exit(-1);
  }
  if(argc>1)
    sscanf(argv[1],"%i",&azimode);
  if(verbose){
    if(azimode)
      fprintf(stderr,"%s: computing azimuth from source to station\n",argv[0]);
    else
      fprintf(stderr,"%s: computing back-azimuth from station to source\n",argv[0]);
  }
  /* 
     convert to radians 
  */
  n=0;
  while(fscanf(stdin,"%lf %lf %lf %lf",&lon0d,&lat0d,&lond,&latd)==4){
    lon0 = lon0d / PIF;
    lat0 = lat0d / PIF;
    colat0d = 90 - lat0d;
    /* 
       sin/cos 
    */
    slat[0]=sin(lat0);clat[0]=cos(lat0);
    
    if(mode == 0){
      /* convert to radian */
      lon = lond/PIF;lat = latd/PIF;
      /* sin/cos */
      slat[1]=sin(lat);clat[1]=cos(lat);
      if(azimode)		/* azimuth  */
				/* back-azimuth */
	azi = atan2(sin(lon0-lon) * clat[0],
		    clat[1] * slat[0] - slat[1] * clat[0] * cos(lon-lon0));
      else
	azi = atan2(sin(lon-lon0) * clat[1],
		    clat[0] * slat[1] - slat[0] * clat[1] * cos(lon0-lon));

      if(azi < 0)
	azi += TWOPI;
      /* distance in radians */
      dist = distance(lon0,lat0,lon,lat,clat[0],clat[1]);
      /* convert to degrees */
      dist *= PIF;azi *= PIF;
    }else{
      colatd=90-latd;
      if(azimode)		/* azimuth */
	stadis(&colatd,&lond,&colat0d,&lon0d,&dist,&azi);
      else			/* back-azimuth */
	stadis(&colat0d,&lon0d,&colatd,&lond,&dist,&azi);
    }
    fprintf(stdout,"%g %g\t%g %g\n",dist,azi,lond,latd);
    n++;
  }
  fprintf(stderr,"%s: processed %i records\n",argv[0],n);
}

/* compute distance on sphere */
double distance(double lon1,double lat1,double lon2,double lat2, 
		double coslat1, double coslat2)
{
  double tmp1,tmp2,tmp3;
  tmp1 = sin((lat1 - lat2)/2.0);
  tmp1 = tmp1 * tmp1;
  
  tmp2 = sin((lon1 - lon2)/2.0);
  tmp2 = tmp2 * tmp2;
  tmp2 *= coslat1;
  tmp2 *= coslat2;

  tmp3 = sqrt(tmp1+tmp2);
  return 2.0*asin(tmp3);
}
 
