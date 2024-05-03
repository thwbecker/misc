#include <stdio.h>
#include <stdlib.h>
#include "bo.h"
#include <math.h>
/* 

program to extract lon lat value ASCII or binary values from binary file
piped into stdin with lon lat z value


usage:


extract_bin_column value [number_of_columns(4)] [select_column_number(3)] ...\
                    [is_double(0)] [ascii_out, 0] [eps, 0.05] [ignore_sign, 0] [print_average, 0]


$Id: extract_bin.c,v 1.4 2005/02/23 03:04:21 becker Exp $

*/
#define OUT_FORMAT "%12g "
int main(int argc,char **argv)
{
  int nin,cn,is_double,i,ascii_out,rbo,n,j,use_abs_eps,ignore_sign,print_average,cnm1;
  float *xflt,rdiff;
  double *xdbl,value,eps,avalue,average,min=1e30,max=-1e30;
  static int verbose = 0;
  ignore_sign = 0;		/* if set, will disregard sign before comparison */
  is_double = 0;		/* is double precision */
  nin = 4;			/* number of fields to read in */
  cn = 3;			/* column to check for value */
  ascii_out =0;			/* ascii output? */
  rbo = 0;			/* reverse byte order */
  
  eps = 5e-02;			/* eps, fractional error . if negative, use absolute */
  
  use_abs_eps = 0;
  
  print_average = 0;
  if(argc < 2){
    fprintf(stderr,"extract_bin_column value [number_of_columns(%i)] [select_column_number(%i)] [is_double(%i)] [ascii_out, %i] [eps, %g] [ignore_sign, %i] [print_average, %i]\n",
	    nin,cn,is_double,ascii_out,eps,ignore_sign,print_average);
    fprintf(stderr,"selects row of binary table with number_of_columns columns, if $(select_column_number) == value\n");
    fprintf(stderr,"eps is the FRACTIONAL error allowed for the match. IF NEGATIVE, WILL USE ABSOLUTE deviations\n");
    exit(-1);
  }
  sscanf(argv[1],"%lf",&value);
  if(argc > 2)
    sscanf(argv[2],"%i",&nin);
  if(argc > 3)
    sscanf(argv[3],"%i",&cn);
  if(argc > 4)
    sscanf(argv[4],"%i",&is_double);
  if(argc > 5)
    sscanf(argv[5],"%i",&ascii_out);
  if(argc > 6)
    sscanf(argv[6],"%lf",&eps);
  if(argc > 7)
    sscanf(argv[7],"%i",&ignore_sign);
  if(argc > 8)
    sscanf(argv[8],"%i",&print_average);

  if(eps < 0){
    use_abs_eps = 1;
    eps = -eps;
  }
  if(verbose){
    fprintf(stderr,"%s: reading %i columns of %s precision binary from stdin\n",
	    argv[0],nin,(is_double)?("double"):("float"));
    fprintf(stderr,"%s: trying to select value %g out of column %i, writing %i remaining fields in %s",
	    argv[0],value,cn,nin-1,(ascii_out)?("ASCII"):("binary"));
    if(ignore_sign)
      fprintf(stderr,", ignoring sign\n");
    else
      fprintf(stderr,"\n");
    if(use_abs_eps)
      fprintf(stderr,"%s: allowing for absolute error of eps %g\n",argv[0],eps);
    else
      fprintf(stderr,"%s: allowing for fractional error of eps %g\n",argv[0],eps);
  }
  if(nin < 1){
    fprintf(stderr,"%s: nin error\n",argv[0]);exit(-1);
  }
  if((cn > nin)||(cn < 1)){
    fprintf(stderr,"%s: cn error\n",argv[0]);exit(-1);
  }

  i=n=0;
  cnm1=cn-1;				
  
  if(ignore_sign)
    value = fabs(value);
  avalue = fabs(value);

  if(!use_abs_eps)
    eps *= avalue;

  average = 0;

  if(is_double){
    xdbl = (double *)malloc(sizeof(double)*(size_t)nin);
    if(!xdbl){fprintf(stderr,"%s: memerror with nin\n",argv[0]);exit(-1);}
    while(fread(xdbl,sizeof(double),nin,stdin)==nin){
      if(rbo)
	flip_byte_order((void *)xdbl,sizeof(double));
      i++;
      if(ignore_sign)
	xdbl[cnm1] = fabs(xdbl[cnm1]);
      rdiff = fabs(xdbl[cnm1] - value);
      if(rdiff < eps){
	n++;
	for(j=0;j < nin;j++){
	  if(j!=cnm1){
	    if(!finite(xdbl[j]))
	      fprintf(stderr,"%s: WARNING: read NaN\n",argv[0]);
	    if(ascii_out)
	      fprintf(stdout,OUT_FORMAT,xdbl[j]);
	    else
	      fwrite((xdbl+j),sizeof(double),1,stdout);
	  }
	}
	average += xdbl[nin-1];
	if(xdbl[nin-1] > max)
	  max = xdbl[nin-1];
	if(xdbl[nin-1] < min)
	  min = xdbl[nin-1];
	if(ascii_out)
	  fprintf(stdout,"\n");
      }
    }
    free(xdbl);
  }else{
    xflt = (float *)malloc(sizeof(float)*(size_t)nin);
    if(!xflt){fprintf(stderr,"%s: memerror with nin\n",argv[0]);exit(-1);}
    while(fread(xflt,sizeof(float),nin,stdin)==nin){
      if(rbo)
	flip_byte_order((void *)xflt,sizeof(float));
      i++;
      if(ignore_sign)
	xflt[cnm1] = fabs(xflt[cnm1]);
      rdiff = fabs(xflt[cnm1] - value);
      if(rdiff < eps){
	n++;
	for(j=0;j < nin;j++){
	  if(j!=cnm1){
	    if(!finite(xflt[j]))
	      fprintf(stderr,"%s: WARNING: read NaN column %i\n",argv[0],j+1);
	    if(ascii_out)
	      fprintf(stdout,OUT_FORMAT,xflt[j]);
	    else
	      fwrite((xflt+j),sizeof(float),1,stdout);
	    
	  }
	}
	average += (double)xflt[nin-1];
	if(xflt[nin-1] > max)
	  max = xflt[nin-1];
	if(xflt[nin-1] < min)
	  min = xflt[nin-1];

	if(ascii_out)
	  fprintf(stdout,"\n");
      }
    }
    free(xflt);
  }
  if(print_average){
    average /= ((double) n);
    if(n)
      fprintf(stderr,"%s: %i x %i in,  %i out for %g(%g) for column %i, average: %g, min: %g, max: %g\n",
	      argv[0],i,nin,n,value,eps,cn,average,min,max);
    else
      fprintf(stderr,"%s: read %i rows, selected none for %g\n",argv[0],i,value);
  }
  if(verbose)
    fprintf(stderr,"%s: read %i rows, selected %i\n",argv[0],i,n);
  return 0;
}
