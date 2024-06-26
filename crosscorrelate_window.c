#include "correlation.h"
#include "correl_nr.h"
/* 

   moving window cross correlation

 */


int main(int argc,char **argv)
{
  COMP_PRECISION *t=NULL,*y1=NULL,*y2=NULL,tcorr,tcorr_max,corr_max,*corr=NULL,dt;
  int i,n,nr,nr2,il,imin,imax,nn;
  COMP_PRECISION window_length = 0.2;	/* window length  */
  int mode = 2;				/* 1: FFT 2: mine 3: mine, Spear */
  COMP_PRECISION rfac = 1;	/* subsample with 1/rfac  */
  
  if(argc < 3){
    fprintf(stderr,"%s file1 file2 [wl, %g] [rfac, %g] [mode, %i]\nreads in two x-y datafiles and computes moving window (fraction wl) cross correlation\n\tresampling at dtmin/rfac\n",
	    argv[0],window_length,rfac,mode);
    fprintf(stderr,"\tif wl < 0, will use actual time interval\n");
    
    exit(-1);
  }
  if(argc>=4)			/* windows length for correlation */
    sscanf(argv[3],"%lf",&window_length);
  if(argc>=5)			/* resampling factor */
    sscanf(argv[4],"%lf",&rfac);
  if(rfac<1){
    fprintf(stderr,"%s: WARNING: sampling larger than unity?! %g\n",argv[0],rfac);
  }
  if(argc >= 6)
    sscanf(argv[5],"%i",&mode);
 
  /* 
     read in two files and interpolate if needed 
  */
  read_two_files_and_interpolate(&t,&y1, &y2,&n,rfac,(argv+1));
  if(window_length < 0){		/* provided as actual time */
    fprintf(stderr,"%s: window length desired in actual time: %g\n",argv[0],-window_length);
    window_length = -window_length/(t[n-1]-t[0]);
  }
  if((window_length < 0)||(window_length>1)){
    fprintf(stderr,"%s: window length error %g (0....1)\n",argv[0],window_length);
    exit(-1);
  }

  nr = (int) ((COMP_PRECISION)n*window_length+0.5);
  if(nr%2!=0)
    nr++;

  dt = t[1]-t[0];
  /* output */
  nr2 = nr/2;
  il = 0;
  imin = il+nr2;imax = n - nr2;
  for(i=imin;i < imax;i++,il++){
    /* time correlation y1_int y2_int */
    tcorr = correlation((y1+i-nr2),(y2+i-nr2),nr,1);

    if(mode == 1){
      compute_correl((y1+i-nr2),(y2+i-nr2),&corr,nr,&nn,0);
    }else{			/* 2 or 3 */
      calc_lag_corr_pedestrian((y1+i-nr2),(y2+i-nr2),&corr,nr,&nn,mode);
    }
    find_max_from_nr_corr(corr,nn,dt,&tcorr_max,&corr_max,mode);
    
    printf("%20.7e\t%20.7f\t%g %g\t\t%20.7e %20.7e\n",
	   t[i],tcorr,tcorr_max,corr_max,y1[i],y2[i]);
  }
  free(t);
  free(y1);
  free(y2);
}

