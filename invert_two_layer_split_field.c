#include "two_layer.h"
/* 



*/
#define RAND ran3

int main(int argc, char **argv)
{
  FILE *fin,*fin2,*foutl,*fouts;
  PREC period = 7;		/* period of SKS  */
  PREC fazi = 0.5;		/* weighting of azimuths compared to
				   delay times from the single station
				   inversion */
  PREC dist_lim = 5;		/* degrees distance */
  PREC dist,adist,min_dist,wsum;
  int i,j,k,min_nsol,max_nsol,min_nb,max_nb;
  char infile[500],outlfile[500],outsfile[500];
  struct sts model;
  struct ind *individual;
  int npop,inv_mode,nelite,imom,idad,nmate;
  PREC felite; 
  struct sl sol;

  model.dw = 0.5;		/* distance weight */
  model.aw = 0.5;	/*	weight of azimuth vs. delay time
				difference */
  model.lw = 0.5;		/* weight of [0] vs. [1] layer */

  model.tlw = 0;		/* top-down variation damping */
  
  model.rseed=-1;

  inv_mode=1;			/* 1: brute force 2: genetic */

  felite=0.1;			/* fraction of solution who can
				   reproduce */

  
  RAND(&model.rseed);			/* init */
  if(argc<2){
    fprintf(stderr,"%s dir [dw, %g] [lw, %g] [aw, %g] [dist_lim, %g] [inv_mode, %i] [tlw, %g] [fazi, %g] [period, %g]\n",
	    argv[0],model.dw,model.lw,model.aw,dist_lim,inv_mode,model.tlw,fazi,period);
    exit(-1);
  }
  sprintf(infile,"%s/stations.dat",argv[1]); /* station info file */
  fin=fopen(infile,"r");
  if(!fin){
    fprintf(stderr,"%s: error with reading station info from %s\n",argv[0],infile);
    exit(-1);
  }
  if(argc>2)
    sscanf(argv[2],FFMT,&model.dw);
  if(argc>3)
    sscanf(argv[3],FFMT,&model.lw);
  if(argc>4)
    sscanf(argv[4],FFMT,&model.aw);
  if(argc>5)
    sscanf(argv[5],FFMT,&dist_lim);
  if(argc>6)
    sscanf(argv[6],"%i",&inv_mode);
  if(argc>7)
    sscanf(argv[7],FFMT,&model.tlw);
  if(argc>8)
    sscanf(argv[8],FFMT,&fazi); /* azimuth weight */
  if(argc>9)
    sscanf(argv[9],FFMT,&period); /* azimuth weight */
  /* 

   */
 
  /* 
     output file name 
  */
  sprintf(outlfile,"inv.d.%4.2f.l.%4.2f.a.%4.2f.tl.%4.2f.dl.%3.1f.fazi.%g.period.%g.log",
	  model.dw,model.lw,model.aw,model.tlw,dist_lim,fazi,period);
  sprintf(outsfile,"sol.d.%4.2f.l.%4.2f.a.%4.2f.tl.%4.2f.dl.%3.1f.fazi.%g.period.%g.dat",
	  model.dw,model.lw,model.aw,model.tlw,dist_lim,fazi,period);

  

  /* 
     stations with distinct locations 
  */
  model.station=(struct st *)malloc(sizeof(struct st));
  model.nstat=0;
  while(fscanf(fin,"%i %s %f %f",&i,
	       model.station[model.nstat].name,
	       &(model.station[model.nstat].lon),&(model.station[model.nstat].lat))==4){
    model.station[model.nstat].sol = (struct sli *)malloc(sizeof(struct sli));
    model.nstat += 1;
    model.station=(struct st *)realloc(model.station,sizeof(struct st)*(model.nstat+1));
  }
  fclose(fin);
  fprintf(stderr,"%s: read %i stations from %s\n",argv[0],model.nstat,infile);

  min_nsol = 1e7;		/* read in the solutions produced by
				   the station-centric two layer
				   fit */
  max_nsol=0;
  for(i=0;i < model.nstat;i++){
    /* solution set */
    sprintf(infile,"%s/fperiod_%g/sols.%i.%g.%g.bin",argv[1],period,i+1,fazi,period);
    fin=fopen(infile,"r");
    if(!fin){
      fprintf(stderr,"%s: error with reading station %s (%i) solution from %s\n",
	      argv[0],model.station[i].name,i+1,infile);
      exit(-1);
    }
    /* basic parameters */
    sprintf(infile,"%s/fperiod_%g/sols.%i.%g.%g.tlfit",argv[1],period,i+1,fazi,period);
    fin2=fopen(infile,"r");
    if(!fin2){
      fprintf(stderr,"%s: error with reading station %s (%i) basic parameters solution from %s\n",
	      argv[0],model.station[i].name,i+1,infile);
      exit(-1);
    }
    if(fscanf(fin2,"%f %f %f %f %f %f %f %f %i %f %f\n",
	      &(model.station[i].sol1l.fazi[0]),&(model.station[i].sol1l.dt[0]),&(model.station[i].sol1l.cost),
	      &(model.station[i].sol2l.cost),
	      &(model.station[i].sol2l.fazi[0]),&(model.station[i].sol2l.dt[0]),
	      &(model.station[i].sol2l.fazi[1]),&(model.station[i].sol2l.dt[1]),
	      &(model.station[i].nobs),&(model.station[i].wfazi),&(model.station[i].period))!=11){
      fprintf(stderr,"%s: basic read error from %s\n",argv[0],infile);
      exit(-1);
    }
    fclose(fin2);
    
    model.station[i].nsol = 0;	/* expect short version of solution,
				   those are sorted by descending
				   misfit */
    while(fread(&(model.station[i].sol[model.station[i].nsol]),
		sizeof(struct sli),1,fin) == 1){
      /* convert reduced chi2 to chi2 by multiplying with degrees of
	 freedom */
      model.station[i].sol[model.station[i].nsol].cost *= (model.station[i].nobs-4);
      //fprintf(stderr,"%s: station %5i solution %6i cost %7.3e\n",
      //argv[0],i,model.station[i].nsol,model.station[i].sol[model.station[i].nsol].cost);
				   
      model.station[i].nsol += 1;
      model.station[i].sol = (struct sli *)realloc(model.station[i].sol,sizeof(struct sli)*(model.station[i].nsol+1));
     }
    fclose(fin);
    model.station[i].sol = (struct sli *)realloc(model.station[i].sol,sizeof(struct sli)*(model.station[i].nsol));
    //fprintf(stderr,"%s: read %i solutions for station %s from %s\n",argv[0],model.station[i].nsol,model.station[i].name,infile);

    if(model.station[i].nsol > max_nsol)
      max_nsol = model.station[i].nsol;
    if(model.station[i].nsol < min_nsol)
      min_nsol = model.station[i].nsol;
  }
  /*  */
  fprintf(stderr,"%s: read solution files, min number: %i, max number: %i computing neighbors for %g deg limit\n",
	  argv[0],min_nsol,max_nsol,dist_lim);

  /* compute spatial relationships */
  min_nb=1e7;max_nb=0;min_dist = 1e10;
  for(i=0;i < model.nstat;i++){
    model.station[i].neighbor = (struct nb *)malloc(sizeof(struct nb));
    if(!model.station[i].neighbor){
      fprintf(stderr,"%s: mem error\n",argv[0]);
      exit(-1);
    }
    model.station[i].nnb=0;
    adist=0;wsum=0;
    for(j=0;j < model.nstat;j++){	/* should be quick even if redundant */
      /*  */
      dist = distance(model.station[i].lon,model.station[i].lat,model.station[j].lon,model.station[j].lat);
      if(i != j){
	if(fabs(dist)<EPS){
	  fprintf(stderr,"%s: ERROR: station %s (%.2f, %2.f) and %s (%.2f, %2.f) have (near) zero distance of %e\n",
		  argv[0],model.station[i].name,
		  model.station[i].lon,model.station[i].lat,model.station[j].name,
		  model.station[j].lon,model.station[j].lat,dist);
	  exit(-1);
	}
	if(dist < dist_lim){	/* found a neighbor */
	  if(dist < min_dist)
	    min_dist = dist;
	  /* store 1/r distance weight */
	  model.station[i].neighbor[model.station[i].nnb].weight = 1/dist;
	  wsum += model.station[i].neighbor[model.station[i].nnb].weight;
	  model.station[i].neighbor[model.station[i].nnb].id = j;
	  model.station[i].nnb += 1;
	  adist += dist;
	  model.station[i].neighbor =
	    (struct nb *)realloc(model.station[i].neighbor,sizeof(struct nb)*(model.station[i].nnb+1));
	  if(!model.station[i].neighbor){
	    fprintf(stderr,"%s: mem error\n",argv[0]);
	    exit(-1);
	  }
	}
      }
    }
    if(model.station[i].nnb == 0){
      fprintf(stderr,"%s: error, station %i %s does not have a neighbor\n",
	      argv[0],i+1,model.station[i].name);
      exit(-1);
    }

    /* normalize weightes */
    for(j=0;j < model.station[i].nnb;j++)
      model.station[i].neighbor[j].weight /= wsum;

    adist /= (PREC)model.station[i].nnb; /* average distance */
    if(model.station[i].nnb > max_nb)
      max_nb = model.station[i].nnb;
    if(model.station[i].nnb < min_nb)
      min_nb = model.station[i].nnb;
    //fprintf(stderr,"%s: found %3i neighbors within %g degress of station %i, avg distance %6.2f deg\n",argv[0],model.station[i].nnb,dist_lim,i+1,adist);
  }
  fprintf(stderr,"%s: found min %i and max %i neighbors, minimum distance: %g\n%s: logfile: %s solfile: %s\n",
	  argv[0],min_nb,max_nb,min_dist,argv[0],outlfile,outsfile);

  /* solutions */
  if(inv_mode == 1){		/* brute force */
    fprintf(stderr,"%s: brute force inversion\n",argv[0]);
    npop = 2;
  }else{			/* genetic */
    npop = 200;			/*  */
    fprintf(stderr,"%s: genetic algorithm with population %i\n",argv[0],npop);
  }
  individual = (struct ind *)malloc(sizeof(struct ind)*npop);
  for(i=0;i<npop;i++){
    individual[i].gene  =  (int *)calloc(model.nstat,sizeof(int));
    if(!individual[i].gene){
      fprintf(stderr,"%s: memerror at pop %i\n",argv[0],npop);
      exit(-1);
    }
  }
  
  /* 
     initialize model 
  */
  foutl = fopen(outlfile,"w");	/* log file */

  /* cost for intial model, station local best solutions (zero
     indices) */
  calc_fitness((individual+0),model);
  switch(inv_mode){
  case 1:			/* brute force random */
    i=0;
    while(i < 1e6){
      perturb_sol_gene(&model,individual[0].gene,individual[1].gene,0.01);
      calc_fitness((individual+1),model);
      if(individual[1].cost  < individual[0].cost){
	individual[0].cost = individual[1].cost;
	individual[0].scost = individual[1].scost;
	individual[0].dcost = individual[1].dcost;
	/* write log of cost functions */
	fprintf(foutl,"%10i t: %10.5e s: %10.5e d: %10.5e\n",i,individual[0].cost,individual[0].scost,individual[0].dcost);fflush(foutl);
	//fprintf(stderr,"%10i t: %10.5e s: %10.5e d: %10.5e\n",i,individual[0].cost,individual[0].scost,individual[0].dcost);
	/* take over better solution */
	memcpy(individual[0].gene,individual[1].gene,model.nstat*sizeof(int));
      }
      i++;
    }
    break;
  case 2:			/* simple genetic algorithm */
    for(i=0;i<npop;i++)	/* initialize with random genes */
      perturb_sol_gene(&model,individual[0].gene,individual[i].gene,1.1); /* all
									     random */
    nelite = (int)npop * felite;
    //nmate = (int) (nelite/2);
    nmate = nelite;
    i=0;
    do{
      /* mate (first iteration will be random) */
      k = nelite;
      while(k < npop){
	imom = (int)(RAND(&model.rseed) * nmate);
	while((idad = (int)(RAND(&model.rseed) * nmate))==imom);;
	//fprintf(stderr,"kid %4i of dad: %4i mom: %4i \n",k-nelite+1,idad+1,imom+1);
	mate(individual[imom],individual[idad],(individual+k),&model,0.01);
	k++;
      }
      for(j=0;j < npop;j++){	/* compute fitness  */
	calc_fitness((individual+j),model);
      }
      /* 
	 sort by fitness 
      */
      qsort(individual,npop, sizeof(struct ind),(int (*)(const void *, const void *))compare_ind);
      /* print progress */
      if(i%5 == 0){
	fprintf(foutl,"%10i t: %10.5e s: %10.5e d: %10.5e\n",i,
		individual[0].cost,individual[0].scost,individual[0].dcost);
	fflush(foutl);
	fprintf(stderr,"%10i t: %10.5e s: %10.5e d: %10.5e\n",
		i,individual[0].cost,individual[0].scost,individual[0].dcost);
      }
      i++;
    }while(i<3500);
    break;
  }
  fclose(foutl);
  /* best solution output */
  fouts = fopen(outsfile,"w");	/* solution file */
  for(i=0;i < model.nstat;i++){
    /* lon lat station name */
    fprintf(fouts,"%g %g %s\t",model.station[i].lon,model.station[i].lat,model.station[i].name);
    /* best overall fit after inversion*/
    sli2sl(model.station[i].sol[individual[0].gene[i]],&sol);
    fprintf(fouts,"%10.5e %5.1f %5.2f %5.1f %5.2f\t",
	    sol.cost,sol.fazi[0],sol.dt[0],sol.fazi[1],sol.dt[1]);
    /* best local fit, without inversion */
    sli2sl(model.station[i].sol[0],&sol);
    fprintf(fouts,"%10.5e %5.1f %5.2f %5.1f %5.2f\n",
	    sol.cost,sol.fazi[0],sol.dt[0],sol.fazi[1],sol.dt[1]);
  }
  fclose(fouts);
  
  return 0;
}
/* 
   mate with mutation rate frac_mod
   
*/
void mate(struct ind dad, struct ind mom,struct ind *kid, struct sts *model,PREC frac_mod)
{
  int i;
  PREC r,f1,f2;
  f2 = 1-frac_mod;
  f1 = f2/2.0;
  
  for(i=0;i < model->nstat;i++){
    r = RAND(&model->rseed);
    if(r < f1){
      kid->gene[i] = dad.gene[i];
    }else if(r < f2){
      kid->gene[i] = mom.gene[i];
    }else{			/* mutation */
      kid->gene[i] = (int)(RAND(&model->rseed) * model->station[i].nsol);
    }
  }
}
  
void calc_fitness(struct ind  *individual, struct sts model)
{
  individual->cost = calc_cost(model,individual->gene,0,
			       &individual->scost,&individual->dcost);
  individual->fitness = 1000./individual->cost;
}
/* 
   frac_mod is the mutation rate 

*/
void perturb_sol_gene(struct sts *model, int *gene,int *genen,PREC frac_mod)
{
  int i;
  for(i=0;i < model->nstat;i++){
    if(RAND(&model->rseed) < frac_mod){
      /* 0.... nsol-1 are possible choices (allowing non
	 modification) */
      genen[i] = (int)(RAND(&model->rseed) * model->station[i].nsol);
    }else{
      genen[i] = gene[i];
    }
  }
}

int compare_ind(struct ind *a, struct ind *b)
{
  if(a->fitness > b->fitness)
    return -1;
  else if(a->fitness == b->fitness)
    return 0;
  else
    return 1;
}

/* compute the total cost of the soli solution model */
PREC calc_cost(struct sts model, int *soli, int verbose,
	       PREC *scost, PREC *dcostt)
{
  PREC dcost_loc[4],dcost[4],weight,dcosttl[2],cost,lvcost[2];
  int i,j,k,inb,itmp;
  struct sli solution,nsolution;
  /* azimuthal scan */
  static int nda=(int)90/TL_DA;
  static int init = 0,damp_two_layer;
  static PREC idw,iaw,ilw,itlw,dw4,tlw4;
  if(!init){
    idw = 1.-model.dw;		/* weight of lateral variation to splitting */
    iaw = 1.-model.aw;		/* weights of azimuth to delay time variations */
    ilw = 1.-model.lw;		/* weight of deep vs. shallow
				   variations */
    dw4 = model.dw/4.0;
    if(model.tlw > 0){		/* damp against variations in each
				   station's top and bottom layer */
      itlw = 1.-model.tlw;
      damp_two_layer = 1;
      tlw4 =  model.tlw/4.0;	/* correct since integer solution
				   difference is scaled by half sigma
				   rather than sigma */
    }else{
      damp_two_layer = 0;
    }
    
    init = 1;
  }
  
  *scost = 0;
  for(j=0;j<4;j++)
    dcost[j] = 0.;
  lvcost[0] = lvcost[1] = 0;	/* */
  for(i=0;i < model.nstat;i++){
    /* this stations current solution */
    solution = model.station[i].sol[soli[i]];
    
    *scost += solution.cost;	/* misfit with splits of this
				   station */

    /* this stations depth variation */
    itmp = abs(solution.iazi[0] - solution.iazi[1]); /* azimuth
							discretized */
    if(itmp > nda)	/* take care of periodicity */
      itmp -= nda;
    lvcost[0] += (PREC)(itmp*itmp); /* chi2 like */
    /* 
       delay time discretized 
    */
    itmp = solution.idt[0] - solution.idt[1];
    lvcost[1] += ((PREC)(itmp*itmp));
    
    for(j=0;j<4;j++)
      dcost_loc[j] = 0.;
    for(j=0;j < model.station[i].nnb;j++){ /* loop through neighbors */
      inb    = model.station[i].neighbor[j].id;
      weight = model.station[i].neighbor[j].weight;
      
      nsolution = model.station[inb].sol[soli[inb]]; /* neighbor solution */
      /* 
	 azimuth mismatch between local and neighbor 
      */
      for(k=0;k<2;k++){					/* bottom and top */
	/* azimuth */
	itmp = abs(solution.iazi[k] - nsolution.iazi[k]);
	if(itmp > nda)	/* take care of periodicity */
	  itmp -= nda;
	dcost_loc[k] += ((PREC)(itmp*itmp))*weight;
	/* delay times */
	itmp = solution.idt[k] - nsolution.idt[k];
	dcost_loc[2+k] += ((PREC)(itmp*itmp))*weight; 
      }
    }
    for(j=0;j<4;j++){
      dcost[j] += dcost_loc[j];
    }
  }
  /* apply weights */
  /* bottom and top weight of delay time vs. azimuth variations */
  dcosttl[0] = iaw*dcost[2]   + model.aw*dcost[0]; /* bottom */
  dcosttl[1] = iaw*dcost[3]   + model.aw*dcost[1]; /* top */
  /* relative weighting of deep vs. shallow variations */
  *dcostt =    ilw*dcosttl[0] + model.lw*dcosttl[1]; /* total distance */
  /* 
     apply a 1/4 correction to the lateral variations
  
  */
  /* 
     splitting (1-dw) vs. lateral variations (dw) 
  */
  cost =       idw*(*scost)   + dw4*(*dcostt);
  if(damp_two_layer){
    /* 
       plus the total two layer variation 
    */
    cost = itlw*cost + tlw4*(iaw*lvcost[1] + model.aw*lvcost[0]);
  }
  if(verbose)
    fprintf(stderr,"tcost: %10.6e scost: %6.3e dcost: %6.3e (%6.3e %6.3e %6.3e %6.3e) tlcost: %6.3e dw: %5.3f aw: %5.3f lw: %5.3f tlw: %5.3f\n",
	    cost,*scost,*dcostt,
	    dcost[0],dcost[1],dcost[2],dcost[3],(iaw*lvcost[1] + model.aw*lvcost[0]),
	    model.dw,model.aw,model.lw,model.tlw);
  
  return cost;
}


/* compute distance in degrees, input in degree */
PREC distance(PREC lon1,PREC lat1,PREC lon2,PREC lat2)
{
  const PREC fac = 57.2957795130823;
  lon1/=fac;
  lat1/=fac;
  lon2/=fac;
  lat2/=fac;
  return distance_rad(lon1,lat1,lon2,lat2,cos(lat1),cos(lat2))*fac;
}


PREC distance_rad(PREC lon1,PREC lat1,PREC lon2,PREC lat2, 
		  PREC coslat1, PREC coslat2)
{
  PREC tmp1,tmp2,tmp3;
  tmp1 = sin((lat1 - lat2)/2.0);
  tmp1 = tmp1 * tmp1;
  
  tmp2 = sin((lon1 - lon2)/2.0);
  tmp2 = tmp2 * tmp2;
  tmp2 *= coslat1;
  tmp2 *= coslat2;

  tmp3 = sqrt(tmp1+tmp2);
  return 2.0*asin(tmp3);
}
 

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-EPS)

PREC ran2(long *idum)
{
  int j;
  long k;
  static int ntabp7 = NTAB + 7;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  PREC temp;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=ntabp7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) 
    *idum += IM1;
  k=idum2/IQ2;
  idum2 = IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) 
    idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) 
    iy += IMM1;
  if ((temp=AM*iy) > RNMX) 
    return RNMX;
  else 
    return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

PREC ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

