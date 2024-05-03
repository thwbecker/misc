#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/* 

read in binary float and turn into a scaled integer, either unsigned or signed 

0: unsigned char      0 ....255
1: signed char     -127 ... 127


$Id: asciifloat2bebin.c,v 1.1 2005/01/07 01:32:14 becker Exp becker $

*/
int main(int argc, char **argv)
{
  int n,type,i,imin,imax;
  unsigned char uchar;
  char schar;
  float *x,xmin,xmax,xrange,xs,oxmin,oxmax;
  static size_t len = sizeof(float);
  if(argc < 2){
    fprintf(stderr,"%s\nread floats as ASCII from stdin and write them in binary integers to stdout\n\n",
	    argv[0]);
    fprintf(stderr,"usage:\n%s type\n\t0: unsigned char 1: signed char 2: symmetric signed char\n",argv[0]);
    exit(-1);
  }
  sscanf(argv[1],"%i",&type);

  n=0;
  x=(float *)malloc(len);
  xmin = 1e20;xmax = -1e20;
  while(fscanf(stdin,"%f",(x+n)) == 1){
    if(x[n] < xmin)xmin = x[n];
    if(x[n] > xmax)xmax = x[n];
    n++;
    x=(float *)realloc(x,(n+1)*len);
    if(!x){
      fprintf(stderr,"%s: memory error during buffering, n=%i\n",
	      argv[0],n);
      exit(-1);
    }
  }


  if(n){
    oxmin = xmin; oxmax = xmax;
    if(type == 2){		/* make symmetric */
      if(fabs(xmin) > fabs(xmax)){
	xmax = -xmin;
      }else{
	xmin = -xmax;
      }
    }
    xrange = xmax - xmin;
    imin = 1e7;imax=-1e7;
    for(i=0;i<n;i++){
      xs = (x[i]-xmin)/xrange;	/* 0..1 */

      if(type == 0){		/* unsigned char */
	uchar = (unsigned char)((int)(xs *255.));
	if((int)uchar < imin)imin = (int)uchar;
	if((int)uchar > imax)imax = (int)uchar;
	fwrite(&uchar,sizeof(unsigned char),1,stdout);
      }else if((type == 1)||(type==2)){	/* signed char */
	schar = (char)((int)(-128. + xs * 255.));
	if((int)schar < imin)imin = (int)schar;
	if((int)schar > imax)imax = (int)schar;
	fwrite(&schar,sizeof(char),1,stdout);
      }
    }
    if(type == 0)
      fprintf(stderr,"%s: %i   unsigned_char data range %8.4e %8.4e scaled range %8.4e %8.4e disc %8.4e irange %i %i\n",
	      argv[0],n,oxmin,oxmax,xmin,xmax,xrange/255,imin,imax);
    else if(type == 1)
      fprintf(stderr,"%s: %i     signed_char data range %8.4e %8.4e scaled range %8.4e %8.4e disc %8.4e irange %i %i\n",
	      argv[0],n,oxmin,oxmax,xmin,xmax,xrange/255,imin,imax);
    else if(type == 2)
      fprintf(stderr,"%s: %i sym signed_char data range %8.4e %8.4e scaled range %8.4e %8.4e disc %8.4e irange %i %i\n",
	      argv[0],n,oxmin,oxmax,xmin,xmax,xrange/255,imin,imax);

  }
  free(x);
  exit(n);
}
