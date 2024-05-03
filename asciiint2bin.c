#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
/* 
   read in n intergers as ASCII and write in binary

*/
int main(int argc, char **argv)
{
  int *i,n;
  if(argc != 1){
    fprintf(stderr,"%s\nread integers as ASCII from stdin and write them in binary to stdout\n\n",
	    argv[0]);
    exit(-1);
  }
  n=0;
  i=(int *)malloc(sizeof(int));
  while(fscanf(stdin,"%i",(i+n)) == 1){
    n++;
    i=(int *)realloc(i,sizeof(int)*(n+1));
    if(!i){
      fprintf(stderr,"%s: memory error during buffering, n=%i\n",
	      argv[0],n);
      exit(-1);
    }

  }
  if(n){
    fwrite(i,sizeof(int),n,stdout);
    fprintf(stderr,"%s: converted %i integers\n",argv[0],n);
  }
  free(i);
  exit(n);
}
