#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//
// read depth in km, produce temp in K for simple half-space cooling
//
// echo depth_in_km | half_space_T age use_adiabat[0/1/2] adiabat_file

double lin_interpolate (double *, double *, int  , double);

int main(int argc, char **argv)
{
  double plate_age,erf_scale,age_s,z,dt,zmin,zmax;
  int use_adiabat,n;
  FILE *in;
  const double Ttop = 273;
  double DT = 1350;	/* difference in C */

  const double sec_per_year = 365.25*24.*60.*60.; // seconds per year
  //kappa = 0.7272e-6;		// diffusivity
  const double kappa = 1e-6;		// diffusivity
  const double adiabat_g = 0.3;		// T/km 100...300 from pyrolite
  double *za,*Ta;
  int na;
  if(argc < 2)
    plate_age=50;
  else
    sscanf(argv[1],"%lf",&plate_age); /* in Myr */
  if(argc < 3)
    use_adiabat=0;
  else
    sscanf(argv[2],"%i",&use_adiabat); /* use an adiabat? */
  if(use_adiabat==1){
    fprintf(stderr,"%s: plate age: %g use_adiabat: %i DT %g agrad: %g, output temp in K\n",
	    argv[0],plate_age,use_adiabat,DT,adiabat_g);
  }else if(use_adiabat==2){
    if(argc >= 4){
      fprintf(stderr,"%s: attempting to read z[km]-T[K] adiabat file from %s\n",argv[0],argv[3]);
      
    }else{
      fprintf(stderr,"%s: attempting to read adiabat but no file given\n",argv[0]);
      exit(-1);
    }
    if(!(in = fopen(argv[3],"r"))){
      fprintf(stderr,"%s: cannot open %s\n",argv[0],argv[3]);
      exit(-1);
    }
    na=0;
    za=(double *)malloc(sizeof(double));
    Ta=(double *)malloc(sizeof(double));
    while(fscanf(in,"%lf %lf",(za+na),(Ta+na))==2){
      na++;
      za=(double *)realloc(za,(na+1)*sizeof(double));
      Ta=(double *)realloc(Ta,(na+1)*sizeof(double));
    }
    fprintf(stderr,"%s: read %i adiabat samples from %s\n",argv[0],na,argv[3]);
    fclose(in);
    DT = lin_interpolate (za,Ta,na,0)-Ttop;
  }else{
    fprintf(stderr,"%s: plate age: %g no adiabat\n",argv[0],plate_age);
  }
  age_s = plate_age * 1.e6 * sec_per_year; // age of plate in s 
  erf_scale = 2.*sqrt(kappa*age_s)/1e3;		// error function scale in km

  //slab_thick = 1.15 * erf_scale * sqrt(2.); // TBL (2.23/2) times sqrt(2) for case D
  zmin=1e20;zmax=-1e20;n=0;
  while(fscanf(stdin,"%lf",&z)==1){
    if(use_adiabat == 1){	/* add simple gradient */
      dt = adiabat_g *z;
      printf("%20.15e\n",Ttop + DT*erf(z/erf_scale) + dt);
    }else if(use_adiabat == 2){	/* from file */
      printf("%20.15e\n",lin_interpolate (za,Ta,na,z)-erfc(z/erf_scale)*DT);
    }else{
      printf("%20.15e\n",Ttop + DT*erf(z/erf_scale));
    }
    if(z<zmin)
      zmin = z;
    if(z>zmax)
      zmax = z;
    n++;
  }
  fprintf(stderr,"%s: processed %i depths from %g to %g km, output in K\n",argv[0],n,zmin,zmax);
  if(use_adiabat==2){
    free(za);free(Ta);
  }
  return 0;
}

double lin_interpolate ( double *x ,double *y, int n , double x1 )
{
  int i,j,n1;
  double fac,fac2;
  double tmp;
  n1 = n-1;
  j=0;
  while( (j < n1) && (x[j] < x1))
    j++;
  if(j == 0)
    j = 1;
  i = j - 1;
  
  fac = (x1 - x[i])/(x[j]-x[i]);
  fac2 = 1.0 - fac;
  tmp=  fac  * y[j] + fac2 * y[i];

  return (tmp);
}
