#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rotate.h"
#include "coordconvention.h"
#include "trig.h"

int read_mat_bis_6(FILE *, COMP_PRECISION [3][3]);
int read_lonlatz(FILE *, COMP_PRECISION [3]);
int print_mat_bis_6(FILE *, COMP_PRECISION [3][3]);

int main(int argc, char **argv)
{
  int n;
  FILE *in1,*in2,*out1;
  COMP_PRECISION r[3],s[3][3],c[3][3];
  if(argc != 4){
    fprintf(stderr,"%s coord_file spherical_matrix_file cart_matrix_file\n",argv[0]);
    fprintf(stderr,"INPUT:\n");
    fprintf(stderr,"reads spherical coordinate system matrix in srr srt srp stt stp spp format, single prec binary from spherical_matrix_file\n");
    fprintf(stderr,"reads spherical coordinate system location in lon(deg) lat(deg) z(km) format, single prec binary from coord_file\n");
    fprintf(stderr,"OUTPUT:\n");
    fprintf(stderr,"prints Cartesian system matrix in sxx sxy sxz syy syz szz format, single prec binary to cart_matrix_file\n" );
    exit(-1);
  }
  in1 = fopen(argv[1],"r");	/* coordinates in lon lat depth format */
  in2 = fopen(argv[2],"r");	/* spherical matrix in six upper
				   triagnle format, single precision
				   binary  */
  out1 = fopen(argv[3],"w");	/* output matrix files */
  if(!in1 || !in2 || !out1 ){
    fprintf(stderr,"%s: could not open %s for locations or %s for matrices or %s for output\n",
	    argv[0],argv[1],argv[2],argv[3]);
    exit(-1);
  }
  fprintf(stderr,"%s: reading lon lat z from %s, six comp spherical tensor single prec binary from %s\n",
	  argv[0],argv[1],argv[2]);
  
  /* loop through matrices */
  n=0;
  while(read_mat_bis_6(in2,s) == 1){ /* matrix */
    if(!read_lonlatz(in1,r)){	     /* location */
      fprintf(stderr,"%s: error: read %i stress tensors, but could not read locations from %s\n",
	      argv[0],n+1,argv[2]);
      exit(-1);
    }
    /* rotate */
    polar_to_cart_mat_at_r3x3(s,c,r);
    /* print cartesian matrix to file */
    print_mat_bis_6(out1,c);
    n++;
  }
  fclose(in1);fclose(in2);fclose(out1);
  fprintf(stderr,"%s: written %i converted matrices to %s\n",argv[0],n,argv[3]);
  return 0;
}

int read_lonlatz(FILE *in, COMP_PRECISION r[3])
{
  int n;
  float tmp[3];
  n = fread(tmp,sizeof(float),3,in);   
  if(n != 3){			/* could not read values */
    return 0;
  }else{
    r[R] = (RADIUS_EARTH_KM-tmp[2])/RADIUS_EARTH_KM;
    r[PHI]=  tmp[0]*PIOVERONEEIGHTY;
    r[THETA]=(90.0-tmp[1])*PIOVERONEEIGHTY;

    return 1;
  }
}

int print_mat_bis_6(FILE *out, COMP_PRECISION s[3][3]){
  int n;
  float tmp[6];
  
  tmp[0] = s[X][X];
  tmp[1] = s[X][Y];
  tmp[2] = s[X][Z];
  tmp[3] = s[Y][Y];
  tmp[4] = s[Y][Z];
  tmp[5] = s[Z][Z];
  n = fwrite(tmp,sizeof(float),6,out);
  if(n != 6){
    return 0;			/*  */
  }else{
    return 1;
  }
}
int read_mat_bis_6(FILE *in, COMP_PRECISION s[3][3]){
  int n;
  float tmp[6];
  n = fread(tmp,sizeof(float),6,in);
  if(n != 6){
    return 0;
  }else{
    s[X][X] = tmp[0];
    s[X][Y] = s[Y][X] = tmp[1];
    s[X][Z] = s[Z][X] = tmp[2];
    s[Y][Y] = tmp[3];
    s[Y][Z] = s[Z][Y] = tmp[4];
    s[Z][Z] = tmp[5];
    return 1;
  }
}
