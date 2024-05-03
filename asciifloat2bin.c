#include <stdio.h>
#include <stdlib.h>
/* 
   read in floats as ASCII and write in binary

*/
int main(int argc, char **argv)
{
  int n;
  float *x;
  if(argc != 1){
    fprintf(stderr,"%s\nread floats as ASCII from stdin and write them in binary to stdout\n\n",
	    argv[0]);
    exit(-1);
  }
  n=0;
  x=(float *)malloc(sizeof(float));
  while(fscanf(stdin,"%f",(x+n)) == 1){
    n++;
    x=(float *)realloc(x,(n+1)*sizeof(float));
    if(!x){
      fprintf(stderr,"%s: memory error during buffering, n=%i\n",
	      argv[0],n);
      exit(-1);
    }
  }
  if(n){
    fwrite(x,sizeof(float),n,stdout);
    fprintf(stderr,"%s: converted %i floats\n",argv[0],n);
  }
  free(x);
  exit(n);
}
