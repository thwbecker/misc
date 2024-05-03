#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "bo.h"
/* 

program to extract single ASCII or binary values from binary file
piped into stdin


usage:


extract_bin_file filename [number_of_columns(1)] [out_column_number(1)] [is_double(0)] [ascii_out, 1]


$Id: extract_bin_file.c,v 1.1 2005/01/06 02:41:45 becker Exp becker $

*/
#define OUT_FORMAT "%12g\n"
int main(int argc,char **argv)
{
  int nin,iout,is_double,i,ascii_out,rbo;
  float *xflt;
  double *xdbl;
  FILE *in;
  is_double = 0;
  nin = 1;
  iout = 1;
  ascii_out =1;
  rbo = 0;			/* reverse byte order */
  if(argc < 2){
    fprintf(stderr,"%s: filename [number_of_columns(%i)] [out_column_number(%i)] [is_double(%i)] [ascii_out(%i)] [reverse_byte_order, %i]\n",
	    argv[0],nin,iout,is_double,ascii_out,rbo);
    exit(-1);
  }
  if(argc > 2)
    sscanf(argv[2],"%i",&nin);
  if(argc > 3)
    sscanf(argv[3],"%i",&iout);
  if(argc > 4)
    sscanf(argv[4],"%i",&is_double);
  if(argc > 5)
    sscanf(argv[5],"%i",&ascii_out);
  if(argc > 6)
    sscanf(argv[5],"%i",&rbo);
  fprintf(stderr,"%s: reading %i columns of %s precision binary from %s, writing column %i to stdout, ",
	  argv[0],nin,(is_double)?("double"):("float"),argv[1],iout);
  if(nin < 1){
    fprintf(stderr,"%s: nin error\n",argv[0]);exit(-1);
  }
  if((iout > nin)||(iout < 1)){
    fprintf(stderr,"%s: iout error\n",argv[0]);exit(-1);
  }
  in = fopen(argv[1],"r");
  if(!in){
    fprintf(stderr,"%s: error opening %s\n",argv[0],argv[1]);
    exit(-1);
  }
  i=0;
  iout--;
  if(is_double){
    xdbl = (double *)malloc(sizeof(double)*nin);
    if(!xdbl){fprintf(stderr,"%s: memerror with nin\n",argv[0]);exit(-1);}
    while(fread(xdbl,sizeof(double),nin,in)==nin){
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
    xflt = (float *)malloc(sizeof(float)*nin);
    if(!xflt){fprintf(stderr,"%s: memerror with nin\n",argv[0]);exit(-1);}
    while(fread(xflt,sizeof(float),nin,in)==nin){
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
  fclose(in);
  return 0;
}

