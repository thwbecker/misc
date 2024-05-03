#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  int n,m,i,j;
  double *x,a,sum;
  FILE *in;
  if(argc!=4){
    fprintf(stderr,"%s m n file\nreads m x n matrix A from stdin and vector x from file\n",
	    argv[0]);
    fprintf(stderr,"writes product A.x to stdout\n");
    exit(-1);
  }
  sscanf(argv[1],"%i",&m);
  sscanf(argv[2],"%i",&n);
  fprintf(stderr,"%s: expecting %i by %i matrix and vector of length %i\n",
	  argv[0],m,n,n);
  if((x=(double *)malloc(sizeof(double)*n))==NULL){
    fprintf(stderr,"%s: memerror with vector\n",argv[0]);
    exit(-1);
  }
  in=fopen(argv[3],"r");
  if(!in){fprintf(stderr,"%s: could not open %s\n",argv[0],argv[3]);exit(-1);}
  for(i=0;i<n;i++){
    if(fscanf(in,"%lf",(x+i))!=1){fprintf(stderr,"%s: read error vector\n",argv[0]);exit(-1);}
    //fprintf(stderr,"%i %g\n",i,x[i]);
  }
  fclose(in);
 
  for(i=0;i<m;i++){
    for(sum=0.0,j=0;j<n;j++)
      if(fscanf(stdin,"%lf",&a)!=1){
	fprintf(stderr,"%s: read error matrix\n",argv[0]);exit(-1);
      }else{
	//fprintf(stderr,"%i %i %g\n",i+1,j+1,a);
	sum += a*x[j];
      }
    fprintf(stdout,"%20.15e\n",sum);
  }
  return 0;
}
