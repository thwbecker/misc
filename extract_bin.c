#include <stdio.h>
#include <stdlib.h>
#include "bo.h"
#include <math.h>
/* 

program to extract single ASCII or binary values from binary file
piped into stdin


usage:


extract_bin [number_of_columns(1)] [out_column_number(1)] [is_double(0)] [ascii_out, 1]


$Id: extract_bin.c,v 1.4 2005/02/23 03:04:21 becker Exp $

*/
#define OUT_FORMAT "%12g\n"
int main(int argc,char **argv)
{
  int nin,iout,is_double,i,ascii_out,rbo;
  float *xflt;
  double *xdbl;
  is_double = 0;
  nin = 1;
  iout = 1;
  ascii_out =1;
  rbo = 0;			/* reverse byte order */
  if(argc > 1)
    sscanf(argv[1],"%i",&nin);
  if(argc > 2)
    sscanf(argv[2],"%i",&iout);
  if(argc > 3)
    sscanf(argv[3],"%i",&is_double);
  if(argc > 4)
    sscanf(argv[4],"%i",&ascii_out);
  if(argc > 5)
    sscanf(argv[5],"%i",&rbo);
  fprintf(stderr,"%s: reading %i columns of %s precision binary from stdin, writing column %i to stdout, ",
	  argv[0],nin,(is_double)?("double"):("float"),iout);
  if(nin < 1){
    fprintf(stderr,"%s: nin error\n",argv[0]);exit(-1);
  }
  if((iout > nin)||(iout < 1)){
    fprintf(stderr,"%s: iout error\n",argv[0]);exit(-1);
  }
  i=0;
  iout--;
  if(is_double){
    xdbl = (double *)malloc(sizeof(double)*(size_t)nin);
    if(!xdbl){fprintf(stderr,"%s: memerror with nin\n",argv[0]);exit(-1);}
    while(fread(xdbl,sizeof(double),nin,stdin)==nin){
      if(rbo)
	flip_byte_order((void *)xdbl,sizeof(double));
      i++;
      if(!finite(xdbl[iout]))
	fprintf(stderr,"%s: WARNING: read NaN\n",argv[0]);
      if(ascii_out)
	fprintf(stdout,OUT_FORMAT,xdbl[iout]);
      else
	fwrite((xdbl+iout),sizeof(double),1,stdout);
    }
    free(xdbl);
  }else{
    xflt = (float *)malloc(sizeof(float)*(size_t)nin);
    if(!xflt){fprintf(stderr,"%s: memerror with nin\n",argv[0]);exit(-1);}
    while(fread(xflt,sizeof(float),nin,stdin)==nin){
      if(rbo)
	flip_byte_order((void *)xflt,sizeof(float));
      i++;
      if(!finite(xflt[iout]))
	fprintf(stderr,"%s: WARNING: read NaN\n",argv[0]);
      if(ascii_out)
	fprintf(stdout,OUT_FORMAT,xflt[iout]);
      else
	fwrite((xflt+iout),sizeof(float),1,stdout);
    }
    free(xflt);
  }
  fprintf(stderr,"read %i rows\n",i);
  return 0;
}
