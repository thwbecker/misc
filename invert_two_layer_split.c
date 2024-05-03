#include "two_layer.h"
/* 
   read in fast_azi dfast_azi dt ddt backazimuth 

   splits and fit a two layer model

   invert_two_layer_split [wfazi, 0.1] [period, 10] [outputfile, itslsol] [single_azimuth, 0] [bot_fazi bot_dt]

   single_azimuth = 0 : fit two layers and have azimuth be variable
   single_azimuth = 1 : fit two layers but have both azimuths be the same 
                        (i.e. only leave dt variable (only for testing purposes))
   single_azimuth = -1: read in bottom layer, fit surface

   uses Kris Walker's implementation of Savage and Silver (1994) via
   Wuestefeld's Matlab implementation of SplitLab, converted to C
  


*/

int main(int argc, char **argv)
{
  PREC mean_fazi,a[2],dt_max,botf_azi,botf_dt;
  struct sl test_sol;
  struct sli *all_soli;		/* short storage */
  PREC period,angle,sin_angle,cos_angle,frac_nsol_dump;
  PREC mean_odt,base_cost,wfazi;
  struct mf chi2,wsum;
  struct sp *obs;
  int nobs,neval,i,nsol_dump,nsol_max,single_azimuth,dof2,dof4,sla_min,sla_max,
    fla_min,fla_max,fld_min,fld_max;
  unsigned short nda,ndt;
  int iazi[2], idt[2];
  FILE *fout;
  char ofname[500],oname[300];
  const PREC pi_fac = 180/M_PI;
  period = 10;			/* 10 s default */
  wfazi = 0.5;			/* weight of azimuth misfit, 0.5 is
				   roughly even */


  single_azimuth = 0;		/* two layers free by default */
  frac_nsol_dump = 0.0025;		/* fraction of solutions to
					   save, set to zero for no
					   output file */
  /* range to scan is hard wired */
  dt_max = 4;
  ndt=(int)dt_max/TL_DT+1;		/* full range, 0...dt_max-dt */
  if(ndt > UCHAR_MAX){
    fprintf(stderr,"%s: ndt (%i) out of range (%i)\n",argv[0],ndt,UCHAR_MAX);
    exit(-1);
  }
  /* azimuthal scan */
  nda=(int)180/TL_DA;
  if(nda > UCHAR_MAX){
    fprintf(stderr,"%s: nda (%i) out of range (%i)\n",argv[0],nda,UCHAR_MAX);
    exit(-1);
  }

  /* arguments */
  if(argc > 1)			/* weights */
    sscanf(argv[1],FFMT,&wfazi);
  if(argc > 2)			/* period */
    sscanf(argv[2],FFMT,&period);
  /*  */
  sprintf (oname,"%s","itslsol"); /* default output filename */
  if(argc > 3)			/* output filename */
    snprintf(oname,200,"%s",argv[3]);
  if(argc > 4)			/* period */
    sscanf(argv[4],"%i",&single_azimuth); /* single fast axis mode */
  switch(single_azimuth){
  case -1:			/* lower layer specified mode */
    if(argc < 6){
      fprintf(stderr,"%s: error, too few parameters\n",argv[0]);
      exit(-1);
    }
    sscanf(argv[5],FFMT,&botf_azi);
    sscanf(argv[6],FFMT,&botf_dt);
    fprintf(stderr,"%s: read in fixed bottom layer parameters azi %g dt %g\n",
	    argv[0],botf_azi,botf_dt);
    break;
  case 0:			/* free */
  case 1:			/* two layers have same azimuth
				   mode */
    ;
    break;
  default:
    fprintf(stderr,"%s: single_azimuth mode %i undefined\n",argv[0],single_azimuth);
    exit(-1);
    ;
  }
  /* 
     init 
  */
  nsol_max = (int)ndt*((int)ndt/2+1)*(int)nda*(int)nda;
  all_soli = (struct sli *)malloc(sizeof(struct sli)*nsol_max);
  if(!all_soli){
    fprintf(stderr,"%s: memory error %i^2 x %i^2\n",argv[0],ndt,nda);
    exit(-1);
  }
  obs = (struct sp *)malloc(sizeof(struct sp));
  
  mean_odt=a[0]=a[1]=0.0;
  wsum.da=wsum.dt=0.0;
  nobs=0;
  /* 
     read in data for station in 

     fazi dfazi dt ddt backazimuth 

     format

  */
  while(fscanf(stdin,F5FMT,&(obs[nobs].phi),&(obs[nobs].wphi),
	       &(obs[nobs].dt),&(obs[nobs].wdt),
	       &(obs[nobs].ba)) == 5){
    /*  */
    if(obs[nobs].wphi < TL_DA)		/* error in phi, mean is 8.32 */
      obs[nobs].wphi = TL_DA;		/* this is usually 5 */
    obs[nobs].wphi = 1/obs[nobs].wphi; /* 1/standard deviation of
					  azimuth */
    /*  */
    if(obs[nobs].wdt < TL_DT)		/* error in dt, mean is 0.25 */
      obs[nobs].wdt = TL_DT;		/* this is usually 0.1 */
    obs[nobs].wdt = 1/obs[nobs].wdt; /* 1/standard deviation of delay
					time */
    /* weight sums */
    wsum.da += obs[nobs].wphi;
    wsum.dt += obs[nobs].wdt;
    /*  */
    angle = obs[nobs].phi*2/pi_fac;
    SCFUNC(angle,&sin_angle,&cos_angle);
    a[0] += sin_angle*obs[nobs].wphi; /* weighted fast axis */
    a[1] += cos_angle*obs[nobs].wphi;
    
    mean_odt += obs[nobs].dt * obs[nobs].wdt; /* weighted delay time */
    /*  */
    nobs++;
    /* make more room */
    obs = (struct sp *)realloc(obs,sizeof(struct sp)*(nobs+1));
  }
  dof2 = nobs - 2;		/* degrees of freedom for two and four
				   layer models */
  dof4 = nobs - 4;
  if(dof4 < 0){
    fprintf(stderr,"%s: error, need at least five observations\n",argv[0]);
    exit(-1);
  }
  
  a[0] /= wsum.da;a[1] /= wsum.da;
  mean_fazi = atan2(a[0],a[1])/2*pi_fac; /* mean fast azimuth */
  mean_odt /= wsum.dt;	/* mean delay time of observations */
  /*  */
  test_sol.fazi[0] = test_sol.fazi[1] = mean_fazi;
  test_sol.dt[0]   = test_sol.dt[1]   = mean_odt/2;
  
  fprintf(stderr,"%s: read %i data, mean fazi: %.1f (da: %.1f) mean dt: %.2f (ddt: %.2f dt_max: %.2f), fazi weight: %g, period: %g\n",
	  argv[0],nobs,mean_fazi,TL_DA,mean_odt,TL_DT,dt_max,wfazi,period);
  
  fprintf(stderr,"%s: da: %g nda %i dt: %g ndt: %i single_azimuth: %i\n",
	  argv[0],TL_DA,(int)nda,TL_DT,(int)ndt,single_azimuth);
  calc_chi2(&chi2,obs,nobs,wsum,test_sol,period); /* chi2 */
  /* weighted reduced chi2 */
  base_cost = cost_function(chi2,wfazi,dof2); /* cost function for
						 single layer given
						 two parameters */
  fprintf(stderr,"%s: single layer reduced chi2: %g (azi), %g (dt) - weighted cost: %g\n",
	  argv[0],chi2.da/((PREC)dof2),chi2.dt/((PREC)dof2),base_cost);
  
  /* 
     parameter space scan 
  */
  neval=0;
  if(single_azimuth == -1){	/* bottom layer is fixed */
    /* make sure within -90...90 range */
    while(botf_azi< -90)
      botf_azi += 180;
    while(botf_azi > 90)
      botf_azi -= 180;
    fla_min = (int)((botf_azi+90)/TL_DA-0.5);
    fla_max = fla_min+1;

    fld_min = (int)(botf_dt/TL_DT-0.5);
    if(fld_min<0)
      fld_min=0;
    fld_max = fld_min+1;
    
    fprintf(stderr,"%s: bottom layer after conversion to increments: azi %g (%i/%i) dt %g (%i/%i)\n",
	    argv[0],-90+(fla_min+0.5)*TL_DA,fla_min,nda,
	    (fld_min+0.5)*TL_DT,fld_min,ndt);

    if((fla_min<0) && (fld_min<0)){
      fprintf(stderr,"%s: out of bounds error at min for phi %i/%i dt %i/%i\n",argv[0],
	      fla_min,0,fld_min,0);

      exit(-1);
    }
    
    if((fla_max > nda) && (fld_max>ndt)){
      fprintf(stderr,"%s: out of bounds error at max for phi %i/%i dt %i/%i\n",argv[0],
	      fla_max,nda,fld_max,ndt);

      exit(-1);
    }
  }else{
    if(single_azimuth == 1)
      fprintf(stderr,"%s: single fast azimuth for both layers\n",argv[0]);
    
    fla_min = 0;		/* bottom azimuth scan */
    fla_max = nda;
    fld_min = 0;		/* bottom layer delay time scan */
    fld_max = ndt;
  }
  for(iazi[0]=fla_min;iazi[0] < fla_max;iazi[0] += 1){ /* bottom layer azimuth */
    if(single_azimuth == 1){
      sla_min=iazi[0];
      sla_max=sla_min+1;
    }else{
      sla_min=0;sla_max=nda;
    }
    for(iazi[1]=sla_min;iazi[1] < sla_max;iazi[1] += 1){ /* top layer azimuth */
      for(idt[0]=fld_min;idt[0] < fld_max;idt[0] += 1){		 /* bottom layer delay time */
	for(idt[1]=0;(idt[0]+idt[1]) < ndt;idt[1] += 1){ /* top layer delay time */
	  for(i=0;i < 2;i++){
	    all_soli[neval].iazi[i] = iazi[i];
	    all_soli[neval].idt[i] =  idt[i];
	  }
	  sli2sl(all_soli[neval],&test_sol);
	  calc_chi2(&chi2,obs,nobs,wsum,test_sol,period); /*  */
	  all_soli[neval].cost = cost_function(chi2,wfazi,dof4);
#ifdef DEBUG
	  fprintf(stderr,"b %5.1f (%3i) %3.1f (%3i) t%5.1f (%3i) %3.1f (%3i) - %8.4f\n",
		  test_sol.fazi[0],iazi[0],
		  test_sol.dt[0],idt[0],
		  test_sol.fazi[1],iazi[1],
		  test_sol.dt[1],idt[1],
		  test_sol.cost);
#endif
	  neval++;
	  if(neval>nsol_max){
	    fprintf(stderr,"%s: error, too many solutions (%i), only provisioned for %i\n",
		    argv[0],neval,nsol_max);
	    exit(-1);
	  }
	}
      }
    }
  }
  all_soli=(struct sli *)realloc(all_soli,neval*sizeof(struct sli));
  fprintf(stderr,"%s: evaluated %i trials (%.1f MB), sorting solutions...\n",
	  argv[0],neval,(PREC)((neval+1)*sizeof(struct sli))/1024/1024);
  /* 
     sort such that best (smallest misfit) is first
  */
  qsort(all_soli,neval, sizeof(struct sli),(int (*)(const void *, const void *))compare_soli);

  /*  */
  fprintf(stderr,"%s: done\n",argv[0]);
  if(frac_nsol_dump > 0){
    /* dump solutions in binary */
    nsol_dump = (int)(frac_nsol_dump*(PREC)neval);
    sprintf(ofname,"%s.%g.%g.bin",oname,wfazi,period);
    fprintf(stderr,"%s: writing %i (%.1f%%) best solutions in binary to %s\n",
	    argv[0],nsol_dump,(PREC)nsol_dump/(PREC)neval*100,ofname);
    fout=fopen(ofname,"w");
    if(!fout){
      fprintf(stderr,"%s: error opening %s\n",argv[0],ofname);
      exit(-1);
    }
    fwrite(all_soli,sizeof(struct sli), nsol_dump,fout);
    fclose(fout);
  }

  sli2sl(all_soli[0],&test_sol);
  /* 

     output is:

    azi_1L dt_1L rchi2_1L rchi2_2L azi_bot dt_bot azi_top dt_top N wfazi period
  */
  fprintf(stdout,"%6.1f %.2f %.4e %.4e %6.1f %.2f %6.1f %.2f %i %g %g\n",
	  mean_fazi,mean_odt,base_cost,
	  test_sol.cost,test_sol.fazi[0],test_sol.dt[0],
	  test_sol.fazi[1],test_sol.dt[1],nobs,wfazi,period);
  free(all_soli);
  return 0;
}



/* this will sort by ascending cost */
int compare_soli(struct sli *a, struct sli *b)
{
  if(a->cost < b->cost)
    return -1;
  else if(a->cost == b->cost)
    return 0;
  else
    return 1;
}




/* 
   cost function is a weighted reduced chi2

   make wfazi between 0 and 1 
*/
PREC cost_function(struct mf chi2,PREC wfazi,int dof)
{
  return (chi2.da*wfazi + chi2.dt*(1-wfazi))/((PREC)dof);
}
/* angular difference in degrees */
PREC cdirdiff(PREC a1, PREC a2, int allow_negative)
{
  PREC da;
   da = a2 - a1;
   if(allow_negative){
     // allow negative (counter vs. clockwise deviations)
     if(da > 0){
       if(da > 180)
	 da -= 180;
       if(da > 90)
	 da -= 180;
     }else{
       if(da<-180)
	 da += 180;
       if(da<-90)
	 da += 180;
     }
   }else{
     // deviations in both ways count the same
     
     if(da < 0)
       da = -da;
     while(da > 180)
       da -= 180;
     if(da > 90)
       da = 180-da;
   }
   return da;
}

/* 

   given a solution and observations compute the fast azimuth and
   delay time misfits, weighted by their uncertainties

   
*/
void calc_chi2(struct mf *chi2,struct sp *obs, int nobs,
	       struct mf wsum,struct sl sol,PREC period)
{
  int i;
  PREC faziapp,dtapp,tmp;
  
  chi2->da = chi2->dt = 0;
  for(i=0;i < nobs;i++){
    //fprintf(stderr,"%g %g %g\n",obs[i].phi,obs[i].dt,obs[i].ba);
    /* compute a two layer model */
    calc_two_layer(sol.fazi,sol.dt,period,obs[i].ba,&faziapp,&dtapp);
    if(finite(faziapp)){	/* chi^2 for azimuth */
      tmp = cdirdiff(faziapp,obs[i].phi,0)*obs[i].wphi;
      chi2->da += tmp*tmp;
    }
    if(finite(dtapp)){		/* chi^2 for dt */
      tmp = fabs(obs[i].dt-dtapp) * obs[i].wdt;
      chi2->dt += tmp*tmp;
    }
  }
}

/* 
   find median of x[n], without messing up x


*/
PREC median(PREC *xo, int n)
{
  int nh,even;
  PREC *x;
  x = (PREC *)malloc(sizeof(PREC)*n);
  if(!x){
    fprintf(stderr,"median: memory error\n");
    exit(-1);
  }
  memcpy(x, xo, sizeof(PREC)*n);
  even = (n%2 == 0)?(1):(0);
  if(even){
    /* even: 1/2(x_{N/2} + x_{N/2+1} */
    nh = n / 2;
    /* call nr_select with offset index */
    return 0.5 * (nr_select(nh,n,(x-1)) + nr_select(nh+1,n,(x-1)));
  }else{
    /* odd: x_{(N+1)/2} */
    return nr_select((n+1)/2,n,(x-1));
  }
  free(x);
}


/* 

numerical recipes routine to find sorted index k of an array with n entries
uses numerical recipes offset index
*/
#define NR_SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

PREC nr_select(int k,int n,PREC *arr)
{
  int i,ir,j,l,mid;
  PREC a,temp;
  
  l=1;
  ir=n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	NR_SWAP(arr[l],arr[ir]);
      }
      return arr[k];
    } else {
      mid=(l+ir) >> 1;
      NR_SWAP(arr[mid],arr[l+1])
	if (arr[l+1] > arr[ir]) {
	  NR_SWAP(arr[l+1],arr[ir]);
	}
      if (arr[l] > arr[ir]) {
	NR_SWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[l]) {
	NR_SWAP(arr[l+1],arr[l]);
      }
      i=l+1;
      j=ir;
      a=arr[l];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	NR_SWAP(arr[i],arr[j]);
      }
      arr[l]=arr[j];
      arr[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}
#undef NR_SWAP

