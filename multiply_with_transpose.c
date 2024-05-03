#include <stdio.h>
#include <stdlib.h>


int main(int argc, char **argv)
{
  int n,m,i,j,k;
  double *a,sum;
  if(argc!=3){
    fprintf(stderr,"%s m n\nreads m x n matrix A from stdin, multiplies it by its transpose\n",
	    argv[0]);
    fprintf(stderr,"and writes the result A^T.A to stdout\n");
    exit(-1);
  }
  sscanf(argv[1],"%i",&m);
  sscanf(argv[2],"%i",&n);
  fprintf(stderr,"%s: expecting %i by %i matrix\n",
	  argv[0],m,n);
  if((a=(double *)malloc(sizeof(double)*n*m))==NULL){
    fprintf(stderr,"%s: memerror with matrix\n",argv[0]);
    exit(-1);
  }
  // read in matrix
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      if(fscanf(stdin,"%lf",(a+i*n+j))!=1){
	fprintf(stderr,"%s: read error matrix\n",argv[0]);exit(-1);
      }
  // calculate transpose
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(sum=0.0,k=0;k<m;k++)
	sum += a[k*n+i]*a[k*n+j];
      fprintf(stdout,"%g ",sum);
    }
    fprintf(stdout,"\n");
  }
      
  
}
