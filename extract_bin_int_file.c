#include <stdio.h>
#include <math.h>
#include "bo.h"
#include <stdlib.h>
/* 

program to extract single ASCII or binary integer values from binary
file piped into stdin


usage:


extract_bin_file filename [number_of_columns(1)] [out_column_number(1)] [ascii_out, 1]


$Id: extract_bin_int_file.c,v 1.2 2005/02/23 03:04:21 becker Exp $

*/
#define OUT_FORMAT "%i\n"
int main(int argc,char **argv)
{
  int nin,iout,i,ascii_out;
  int *xint;
  FILE *in;
  nin = 1;
  iout = 1;
  ascii_out =1;
  if(argc < 1){
    fprintf(stderr,"%s: filename [number_of_columns(%i)] [out_column_number(%i)] [ascii_out(%i)]\n",
	    argv[0],nin,iout,ascii_out);
    exit(-1);
  }
  
  if(argc > 1)
    sscanf(argv[1],"%i",&nin);
  if(argc > 2)
    sscanf(argv[2],"%i",&iout);
  if(argc > 3)
    sscanf(argv[4],"%i",&ascii_out);
  fprintf(stderr,"%s: reading %i columns of integer binary from stdin, writing column %i to stdout, ",
	  argv[0],nin,iout);
  if(nin < 1){
    fprintf(stderr,"%s: nin error\n",argv[0]);exit(-1);
  }
  if((iout > nin)||(iout < 1)){
    fprintf(stderr,"%s: iout error\n",argv[0]);exit(-1);
  }

  in = fopen(argv[1],"r");
  if(!in){
    fprintf(stderr,"%s: error opening %s\n",
	    argv[0],argv[1]);
    exit(-1);
  }

  i=0;
  iout--;

  xint = (int *)malloc(sizeof(int)*nin);
  if(!xint){fprintf(stderr,"%s: memerror with nin\n",argv[0]);exit(-1);}
  while(fread(xint,sizeof(int),nin,in)==nin){
    i++;
    if(!finite(xint[iout]))
      fprintf(stderr,"%s: WARNING: read NaN\n",argv[0]);
    if(ascii_out)
      fprintf(stdout,OUT_FORMAT,xint[iout]);
    else
      fwrite((xint+iout),sizeof(int),1,stdout);
  }
  free(xint);
  free(in);
  fprintf(stderr,"read %i rows\n",i);
  return 0;
}
