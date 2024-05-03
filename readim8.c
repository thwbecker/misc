#include <stdlib.h>
#include <stdio.h>
#define TYPE unsigned char

void main(int argc,char **argv)
{
  int n;
  TYPE in;
  fprintf(stderr,"%s: reading %i byte values from stdin\n",argv[0],sizeof(TYPE));

  n=0;
  while(fread(&in,sizeof(TYPE),1,stdin)==1){
    printf("%i\n",(int)in);
    n++;
  }
  fprintf(stderr,"%s: read %i values\n",argv[0],n);
}
