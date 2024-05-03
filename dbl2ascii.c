#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv)
{
  int n,i=0,j;
  double *x;
  if(argc!=2){
    fprintf(stderr,"%s number_of_doubles_in_row\n",argv[0]);
    fprintf(stderr,"\tconverts rows of binary data to ascii\n"); 
    fprintf(stderr,"\treads from stdin, writes to stdout\n");
    exit(-1);
  }
  sscanf(argv[1],"%i",&n);
  fprintf(stderr,"%s: expecting %i items in each row\n",
	  argv[0],n);
  x=(double *)malloc(sizeof(double)*n);
  while(fread(x, sizeof(double),n,stdin)==n){
    for(j=0;j<n;j++)
      fprintf(stdout,"%g ",x[j]);
    fprintf(stdout,"\n");
    i++;
  }

  fprintf(stderr,"%s: %i lines read\n",argv[0],i);
  free(x);
  return 0;
}
