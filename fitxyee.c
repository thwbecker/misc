#include "fitxyee.h"


/*

  fit a line through x y data with uncertainties in both x and y

  $Id: fitxyee.c,v 1.2 2005/07/04 16:47:43 becker Exp becker $

  From Numerical Recipes in C, p.668


  major modifications:

  - using data structure instead of x,y,sigx,sigy
  - removed global variables and put those into fit structure

*/


int main(int argc, char **argv)
{
  int nrp,verbose=1,mode,i,twolinear,i2left,i2best,ilim,
    offset,skipped;
  PRECISION a[2],b[2],siga[2],sigb[2],chi2[2],
    q[2],dum1,varx,vary,xeps,yeps,mx,my,minchi2,chi2t;
  int wfile;
  FILE *out;
  struct dp *data;
  mode = 0;
  twolinear = 0;
  if((argc > 4)||((argc==2)&&(strcmp(argv[1],"-h")==0))){
    fprintf(stderr,"%s [mode, 0] [twolinear, 0] [fit_function]\n\tfits a straight line to data with errors in x and y\n",
	    argv[0]);
    fprintf(stderr,"\treads x y sigx sigy from stdin\n");
    fprintf(stderr,"\tnumerics from numerical recipes, chapter 15.3\n");
    fprintf(stderr,"\toutput is a b siga sigb chi2 q\n");
    fprintf(stderr,"\tif mode is set to 1, will assume both x and y have fixed sigmas\n");
    fprintf(stderr,"\tif mode is set to 2, will assume both x and y have relative sigmas\n");
    fprintf(stderr,"\tfor both mode 1 and 2, will read x y couples only\n\n");
    fprintf(stderr,"\tif twolinear is set to unity, will find two best-fit linear fits\n");
    fprintf(stderr,"\tif a third argument is given, will print a fitting function to that file\n\n");
    exit(-1);
  }
  if(argc > 1){// select mode

    sscanf(argv[1],"%i",&mode);
  }
  if(argc > 2){// select mode
    sscanf(argv[2],"%i",&twolinear);
  }
  if(argc > 3){			/* open file */
    wfile = 1;
    out = fopen(argv[3],"w");
    if(!out){
      fprintf(stderr,"%s: couldn't open %s\n",argv[0],argv[3]);
      exit(-1);
    }
    fprintf(stderr,"%s: writing to %s\n",
	    argv[0],argv[3]);
  }else{
    wfile=0;
  }

  data=(struct dp *)malloc(sizeof(struct dp));
  nrp=0;skipped=0;
  if(mode == 0){// read in sigma
    if(verbose)
      fprintf(stderr,"%s: reading x y sig_x sig_y\n",argv[0]);
    while(fscanf(stdin,FOURFORMAT,&(data+nrp)->x,&(data+nrp)->y,
		 &(data+nrp)->sigx,&(data+nrp)->sigy) == 4){
      if(finite((data+nrp)->x)&&finite((data+nrp)->y)&&finite((data+nrp)->sigx)&&finite((data+nrp)->sigy)){
	nrp++;
	data=(struct dp *)realloc(data,sizeof(struct dp)*(nrp+1));
	if(!data){
	  fprintf(stderr,"%s: memerror, n: %i\n",argv[0],nrp);
	  exit(-1);
	}
      }else{
	skipped++;
      }
    }
  }else{// no sigmas given
    if(verbose)
      fprintf(stderr,"%s: reading x y, assuming sig_x=sig_y=1.0\n",
	      argv[0]);
    while(fscanf(stdin,TWOFORMAT,&(data+nrp)->x,
		 &(data+nrp)->y) == 2){
      if(finite((data+nrp)->x)&&finite((data+nrp)->y)){
	/* assign constant uncertainty */
	data[nrp].sigx = data[nrp].sigy = 1.0;
	nrp++;
	data=(struct dp *)realloc(data,sizeof(struct dp)*(nrp+1));
	if(!data){
	  fprintf(stderr,"%s: memerror, n: %i\n",argv[0],nrp);
	  exit(-1);
	}
      }else{
	skipped++;
      }
    }
  }
  if(mode == 2){
    /* 
       
    scale uncertainties
    
    */
    mx=my=0.0;
    for(i=0;i<nrp;i++){		/* find mean */
      mx += data[i].x;my+=data[i].x;
    }
    mx/=nrp;my/=nrp;
    varx=vary=0.0;		/* find var from mean */
    for(i=0;i<nrp;i++){
      dum1 = data[i].x-mx;dum1*=dum1;varx+=dum1;
      dum1 = data[i].y-my;dum1*=dum1;vary+=dum1;
    }
    varx = sqrt(varx/nrp);    vary = sqrt(vary/nrp);
    xeps=varx*0.1;
    yeps=vary*0.1;
    /* scale uncertainties with relative size */
    for(i=0;i<nrp;i++){
      dum1 = fabs(data[i].x) * 0.05;
      if(dum1 < xeps)data[i].sigx = xeps;else data[i].sigx=dum1;
      dum1 = fabs(data[i].y) * 0.05;
      if(dum1 < yeps)data[i].sigy = yeps;else data[i].sigy=dum1;
    }
  }
  /* sort data */
  qsort(data,nrp,sizeof(struct dp),
	(int (*)(const void *, const void *))comparef);
  /*  */
  if(verbose)
    fprintf(stderr,"%s: read n: %i data points, xmin: %g xmax: %g, %i were not finite/unusable\n",
	    argv[0],nrp,data[0].x,data[nrp-1].x,skipped);
  if(!nrp){
    fprintf(stderr,"%s: Exiting, no data.\n",argv[0]);exit(-1);}
  if(!twolinear){
    //
    // call fitting routine
    //
    fitexy((data-1),nrp,a,b,siga,sigb,chi2,q);
    fprintf(stderr,"%s: total chi2: %14.7e nchi2: %14.7e\n",
	    argv[0],chi2[0],sqrt(chi2[0]/(nrp-2)));
    fprintf(stderr,"%s: writing a b siga sigb chi2 q n chi2_hat\n",argv[0]);
    fprintf(stdout,"%14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %8i %14.7e\n",
	    a[0],b[0],siga[0],sigb[0],chi2[0],q[0],nrp,sqrt(chi2[0]/(nrp-2)));
    if(wfile){
      for(i=0;i<nrp;i++)
	fprintf(out,"%11g %11g\n",
		data[i].x,a[0] + data[i].x * b[0]);
      fclose(out);
    }
  }else{
    /* 
       find two, bilinear best-fit lines 
    */
    minchi2 = 1e20;
    if(nrp < 6){
      fprintf(stderr,"%s: error, need at least six points for bilinear fit\n",
	      argv[0]);
      exit(-1);
    }
    /* scan through ...3 ... nrp-3... */
    ilim=nrp-3;
    for(i2left=3;i2left<ilim;i2left++){	/* scan through all */
      /* left fit */
      offset = -1;
      fitexy((data+offset),i2left,a,b,siga,sigb,chi2,q);
      /* right fit */
      offset = -1+i2left;
      fitexy((data+offset),nrp-i2left,(a+1),(b+1),(siga+1),
	     (sigb+1),(chi2+1),(q+1));
      /* 
	 total normalized misfit is used to decide best 
	 division
      */
      chi2t = chi2[0]/(i2left-2) + chi2[1]/(nrp-i2left-2);
      if(chi2t < minchi2){
	minchi2 = chi2t;
	i2best = i2left;
      }
    }
    /* compute best fits */
    offset = -1;
    fitexy((data+offset),i2best,a,b,siga,sigb,chi2,q);
    offset = -1+i2best;
    fitexy((data+offset),nrp-i2best,(a+1),(b+1),(siga+1),
	   (sigb+1),(chi2+1),(q+1));
    /* misfits */
    chi2t = chi2[0] + chi2[1];
    fprintf(stderr,"%s: total chi2: %14.7e nchi2: %14.7e\n",
	    argv[0],chi2t,sqrt(chi2t/(nrp-4)));
    fprintf(stderr,"%s: writing a_i b_i siga_i sigb_i chi2_i q_i n_i x_min x_max\n",
	    argv[0]);
    fprintf(stdout,"%14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %8i %11g %11g\n",
	    a[0],b[0],siga[0],sigb[0],chi2[0],q[0],i2best,data[0].x,data[i2best-1].x);
    fprintf(stdout,"%14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %8i %11g %11g\n",
	    a[1],b[1],siga[1],sigb[1],chi2[1],q[1],nrp-i2best,data[i2best].x,data[nrp-1].x);
    if(wfile){
      for(i=0;i<i2best;i++)
	fprintf(out,"%11g %11g\n",
		data[i].x,a[0] + data[i].x * b[0]);
      for(i=i2best;i<nrp;i++)
	fprintf(out,"%11g %11g\n",
		data[i].x,a[1] + data[i].x * b[1]);
      
      fclose(out);
    }  
  }
  return 0;
}
