#include <stdlib.h>

#include <stdio.h>
/* 
   read in n unsinged intergers as ASCII and write in binary

*/
int main(int argc, char **argv)
{
  unsigned int i,n;
  if(argc != 1){
    fprintf(stderr,"%s\nread integers as ASCII from stdin and write them in binary to stdout\n\n",
	    argv[0]);
    exit(-1);
  }
  n=0;
  while(fscanf(stdin,"%i",&i)==1){
    fwrite(&i,sizeof(unsigned int),1,stdout);
    n++;
  }
  fprintf(stderr,"%s: converted %i unsinged integers\n",argv[0],n);
  exit(n);
}
