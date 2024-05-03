#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sacio.h>
/* 

   simple program to read SAC files


 */
#define MAX 100000

#define CHECK_ERROR {if(nerr > 0){fprintf(stderr,"%s: error %i\n",argv[0],nerr);exit(nerr);}}
int main(int argc, char **argv)
{
  float yarray[ MAX ] , beg , del ,delta,slon,slat,elon,elat,edist;
  int nlen , nerr , max = MAX ,mode=1,i;
  
  if(argc < 2){
    fprintf(stderr,"%s file [mode, %i]\n",argv[0],mode);
    fprintf(stderr,"mode=1: x-y\n");
    fprintf(stderr,"mode=2: slon slat elon elat delta\n");
    fprintf(stderr,"mode=3: great circle distance\n");
    exit(-1);
  }
  if(argc>2)
    sscanf(argv[2],"%i",&mode);
  
  rsac1(argv[1], yarray, &nlen, &beg, &del, &max, &nerr, strlen(argv[1]) ) ;CHECK_ERROR; 
  getfhv ( "DELTA" , &delta , &nerr , 5 ) ; CHECK_ERROR; /* get spacing */
  getfhv ( "STLO" , &slon , &nerr , 5 ) ; CHECK_ERROR; 
  getfhv ( "STLA" , &slat , &nerr , 5 ) ; CHECK_ERROR; 
  getfhv ( "EVLO" , &elon , &nerr , 5 ) ; CHECK_ERROR; 
  getfhv ( "EVLA" , &elat , &nerr , 5 ) ; CHECK_ERROR; 
  getfhv ( "GCARC", &edist, &nerr, 5); CHECK_ERROR;
  switch(mode){
  case 1:			/* x - y */
    for(i=0;i<nlen;i++)
      fprintf(stdout,"%15.7e %15.7e\n",i*delta,yarray[i]);
    break;
  case 2:			/* parameters */
    fprintf(stdout,"%g %g %g %g %g\n",slon,slat,elon,elat,delta);
    break;
  case 3:			/* distance */
    fprintf(stdout,"%15.7e\n",edist);
    break;
  default:
    fprintf(stderr,"%s: mode %i undefined\n",argv[0],mode);
    exit(-1);
  }

}
