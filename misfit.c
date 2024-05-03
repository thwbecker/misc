#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc,char **argv)
{
  int n,m,i,j;
  float *a,*b,*x,misfit,norm,sum;
  if(argc!=3){
    fprintf(stderr,"%s m n\n",argv[0]);
    fprintf(stderr,"reads a.dat (mxn) x.dat(n) b.dat(m) from stdin\n");
    fprintf(stderr,"prints misfit and model norm\n");
    exit(-1);
  }
  sscanf(argv[1],"%i",&m);
  sscanf(argv[2],"%i",&n);
  a=(float *)malloc(sizeof(float)*m*n);
  b=(float *)malloc(sizeof(float)*m);
  x=(float *)malloc(sizeof(float)*n);
  if(!x || !b || !a){
    fprintf(stderr,"memerror\n");
    exit(-1);
  }
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      if(fscanf(stdin,"%f",(a+i*n+j))!=1){
	fprintf(stderr,"read error a\n");
	exit(-1);
      }
  for(j=0;j<n;j++)
    if(fscanf(stdin,"%f",(x+j))!=1){
      fprintf(stderr,"read error x\n");
      exit(-1);
    }
  for(i=0;i<m;i++)
    if(fscanf(stdin,"%f",(b+i))!=1){
      fprintf(stderr,"read error b\n");
      exit(-1);
    }
  for(misfit=0.0,i=0;i<m;i++){
    for(sum=0.0,j=0;j<n;j++)
      sum += a[i*n+j]*x[j];
    sum -= b[i];
    misfit += sum*sum;
  }
  misfit = sqrt(misfit);

  for(norm=0.0,j=0;j<n;j++)
    norm += x[j]*x[j];
  norm=sqrt(norm);

  printf("%g %g\n",misfit,norm);

  free(a);free(b);free(x);
  return 0;
}
