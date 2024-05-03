/* 

using GMT grd files in 3-D for tomographic models, detect
plumes 

$Id: plume_detect.c,v 1.2 2006/12/09 22:03:18 becker Exp becker $

*/
#include "ggrd_grdtrack_util.h"
#include <limits.h>
#include <stdlib.h>
#include <malloc.h>

#define NXY(i,j) ((i)*nlon + (j))

/* use unsigned char */
#define PCTYPE unsigned char
#define NMAX_PLUME 127
//#define NMAX_PLUME UCHAR_MAX


#define R_EARTH 6371.0

/* pixcel, with two integer coord */
struct pxl{
  unsigned int c[2];		/* j (0...nlat-1), i (0...nlon-1)  */
};
/* 

   plume structure 

*/
struct plm{
  int nz,*iz,			/* number of depth levels, and depth
				   level code */
    *npixel;			/* number of pixels per layer */
  struct pxl **pixels;		/* pixels[nz] */
  int *pnr,**pname;		/* patch numbers[nz] and
				   names[nz][pnr[nz]] */
  int mean_npixel;
  GGRD_CPREC *area,mean_area,total_area,zmin,zmax,volume;
};
/* 
   layer structure  
*/
struct lyr{
  GGRD_CPREC z,pa;			/* depth, total patch area*/
  GGRD_CPREC tmean,trms,tmax;	/* mean, rms, and max */
  PCTYPE *patch_code;	/* all codes, which ones assigned */
  int npatch;
  /* per patch */
  unsigned int *ppixcount;
  struct pxl *pcentroid;
  GGRD_CPREC *parea;
};

void merge_plumes(struct plm **, int *,
		  int , int ,struct lyr *, int , int );

int check_for_patch(int , struct plm *,
		    int ,int , int, struct lyr *, 
		    int , int ,GGRD_CPREC ,GGRD_CPREC *, 
		    GGRD_CPREC *,GGRD_CPREC );
void lonlat2xyz(GGRD_CPREC , GGRD_CPREC ,GGRD_CPREC *);
void xyz2lonlat(GGRD_CPREC *,GGRD_CPREC *, GGRD_CPREC *);
void make_new_plume(struct plm **,int *);
void init_plume(struct plm *);
void clear_plume(struct plm *);
void free_plume(struct plm **);
void add_patch_to_plume(int , int , int,struct lyr *, struct plm *, int ,
			int );
GGRD_CPREC int_dist(int ,int ,int ,int ,GGRD_CPREC *,GGRD_CPREC *,
	       GGRD_CPREC );
void compute_mean_plume_properties(struct plm *, int,GGRD_CPREC *, GGRD_CPREC ,struct lyr *,
				   GGRD_CPREC);
void print_plume(struct plm * ,FILE *,GGRD_CPREC *, GGRD_CPREC *,struct lyr *);
int plm_compare(struct plm *, struct plm *);
void add_plume_to_px_assign(struct plm *,int , int , struct lyr *,
			    PCTYPE *);
GGRD_CPREC nd_r_of_zkm(GGRD_CPREC );
GGRD_CPREC compute_total_pix_assigned_vol(int ,int ,int ,struct lyr *,
					  PCTYPE *,GGRD_CPREC *, GGRD_CPREC ,
					  GGRD_CPREC );
void print_patch(struct lyr *,FILE *,GGRD_CPREC *, GGRD_CPREC *,int , int );
PCTYPE nbr_code(int, int, int, int, struct lyr *,int);

int main(int argc, char **argv)
{
  GGRD_CPREC *lon,*lat,ws,*w,frac,*temp,dlon,dlat,dt,tcut,tvol,
    dz,dzw,*xp=NULL,xloc[3],tmplon,tmplat,rmax,marea,dfac,min_length,
    max_depth,scale,tmp,zrange,r2,r,rvol;
  struct ggrd_gt t3d[1];
  double value;
  static double pif = GGRD_PI/180;	
  int i,j,k,nlon,nlat,os,iz,nn,ii,jj,minlo,maxlo,minla,maxla,found,
    midlo,midla,istart,iloop,jstart,jstop,jloop,nz,apatch,nplume,nzero,
    nplume_select,old_found,nfound,loop,iiloop;
  struct lyr *layer;
  struct plm *plume = NULL;
  char filename[300];
  FILE *out,*out2;
  PCTYPE *passigned,pfound,plow,phigh;
  static GGRD_CPREC zbot = 2750, ztop = 300;

#ifndef PD_DOWN			/* plume_detect mode */
  static int upward = 1;	/* 
				   
				0: search top down
				1: search bottom up

				*/
#else  /* plume_detect_d */
  static int upward = 0;
#endif
  /* if this is set to zero, will only allow upward branching for the
     upward method
  */
  static int allow_downward_branch = 0;

  /* 
     plume detection settings 
  */
  frac = .5;			/* cutoff level  */

  dz = 50;			/* depth spacing */

  dfac = 3;			/* how many depth spacings to allow
				   vertically? */
  min_length = 500;		/* minimum length for output */
  max_depth = zbot;		/* only output if shallowest depth is smaller than */

  scale = -1;			/* scale for grids 	T = scale * value */

  /* parameterization settings */
  dlon = dlat = 1;		/* spacing of interpolation grid */

  if(argc < 3){
    fprintf(stderr,"usage:\n\n%s prefix(grd) depth_file [frac, %g] [dfac, %g] [min_length, %g] [max_depth, %g] [scale, %g]\n",
	    argv[0],frac,dfac,min_length,max_depth,scale);
    
    exit(-1);
  }

  if(argc > 3)
    sscanf(argv[3],"%lf",&frac);
  if(argc > 4)			/*  */
    sscanf(argv[4],"%lf",&dfac);
  if(argc > 5)			/*  */
    sscanf(argv[5],"%lf",&min_length);
  if(argc > 6)			/*  */
    sscanf(argv[6],"%lf",&max_depth);
  if(argc > 7)			/*  */
    sscanf(argv[7],"%lf",&scale);

  if((max_depth < ztop)||(max_depth>zbot)){
    fprintf(stderr,"%s: error: max_depth (%g) should be between ztop: %g and zbot: %g\n",
	    argv[0],max_depth,ztop,zbot);
    exit(-1);
  }
  fprintf(stderr,"%s: prefix: %s dfile: %s, using frac: %g dfac: %g at dz: %g\n\t\tmin_length: %g max_depth: %g scale: %g\n",
	  argv[0],argv[1],argv[2],frac,dfac,dz,min_length,max_depth,scale);
  if(frac < 0)
    fprintf(stderr,"%s: using mean + %g * RMS for cutoff\n",argv[0],-frac);
  else
    fprintf(stderr,"%s: using mean + %g * (tmax-tmean) for cutoff\n",argv[0],frac);
  if(upward){
    fprintf(stderr,"\n%s: using old upward detector algorithm, %s downward branching\n\n",
	    argv[0],(allow_downward_branch)?("with"):("without"));
  }else{
    fprintf(stderr,"\n%s: using new downward detector algorithm with double counting\n\n",
	    argv[0]);
  }
  rmax = dfac * dz;			/* max distance, in km */
  /* number of samples */
  nlon = 360/dlon;
  nlat = 180/dlat;
  nz = (zbot-ztop)/dz;
  nn = nlon * nlat;

  /* allocate arrays */
  temp = (GGRD_CPREC *)malloc(sizeof(GGRD_CPREC)*nn); /* temperatures */
  w = (GGRD_CPREC *)malloc(sizeof(GGRD_CPREC)*nlat); /* weights(lat) */
  lon = (GGRD_CPREC *)malloc(sizeof(GGRD_CPREC)*nlon); /* lon coord */
  lat = (GGRD_CPREC *)malloc(sizeof(GGRD_CPREC)*nlat); /* lat coord */
  layer = (struct lyr *)calloc(nz,sizeof(struct lyr));
  if(!lon || !lat || !temp || !w || !layer){
    fprintf(stderr,"%s: mem error\n",argv[0]);exit(-1);}

  /* 
     read in GMT grids 
  */
  t3d[0].init = FALSE;
  if(ggrd_grdtrack_init_general(TRUE,argv[1], argv[2],"-fg",t3d,FALSE,FALSE,FALSE)){
    fprintf(stderr,"%s: ggrd init error\n",argv[0]);
    exit(-1);
  }
  /* 
     loop through layers bottom up
  */
  for(iz=0;iz < nz;iz++){
    /* 

    main depth loop

    */

    if(iz==0)
      layer[iz].z = zbot;
    else
      layer[iz].z = layer[iz-1].z - dz;
    /* 
       allocate memory for this layer and set to zero
    */
    layer[iz].patch_code = (PCTYPE *)calloc((size_t)nn,sizeof(PCTYPE));
    if(!layer[iz].patch_code){
      fprintf(stderr,"%s: mem error patch code \n",argv[0]);
      exit(-1);
    }
    layer[iz].tmax = -1e20;
    /*
       layer 
    */
    if(iz == 0)
      ws = 0.0;			/* init weight sum */
    /* 

    read in temps and determine mean as well as max

    */
    for(i=os=0;i < nlat;i++){ /* lat loop */
      if(iz == 0){		/* init */
	if(i==0)
	  lat[i] = -90+dlat/2.0;
	else
	  lat[i] = lat[i-1]+dlat; /* coords */
	w[i] = cos(pif * lat[i]); /* weight for this lat */
      }
      for(j=0;j < nlon;j++,os++){ /* lon loop */
	if(iz==0){		/* init */
	  ws += w[i];
	  if(j==0)
	    lon[j] = 0.0;	/*  */
	  else
	    lon[j] = lon[j-1] + dlon;
	}
	/* determine velocity */
	if(!ggrd_grdtrack_interpolate_xyz(lon[j],lat[i],layer[iz].z,t3d,&value,TRUE)){
	  fprintf(stderr,"%s: interpolate error: %g %g %g\n",argv[0],
		  lon[j],lat[i],layer[iz].z);
	  exit(-1);
	}
	/* go from velocity to "temperature" */
	temp[os] = scale * value;

	/* layer mean */
	layer[iz].tmean += w[i] * temp[os];
	/* layer max */
	if(layer[iz].tmax < temp[os])
	  layer[iz].tmax = temp[os];
      }
    }
    layer[iz].tmean /= ws;
    /* 
       compute RMS
    */
    layer[iz].trms = 0.0;
    for(i=os=0;i<nlat;i++)
      for(j=0;j<nlon;j++,os++){
	tmp = temp[os] - layer[iz].tmean;
	layer[iz].trms +=  w[i] * tmp*tmp;
      }
    layer[iz].trms = sqrt(layer[iz].trms /ws);
  
     /* 
       cutoff temperatures 
    */
    if(frac > 0){
      /* use fraction of range from mean to max */
      dt = layer[iz].tmax - layer[iz].tmean;
      tcut = frac*dt + layer[iz].tmean;	/* cutoff temperature */
    }else{
      /* use multiples of the RMS  */
      tcut = layer[iz].tmean - frac * layer[iz].trms;
      if(tcut > layer[iz].tmax){
	fprintf(stderr,"%s: adjusting: %g times the RMS is outside the maximum of %g\n",
		argv[0],-frac,layer[iz].tmax);
	tcut = layer[iz].tmax;
      }
    }
    
    //fprintf(stderr,"%s: z: %g tmean: %g tmax: %g rms: %g\n",argv[0],layer[iz].z,layer[iz].tmean,layer[iz].tmax,layer[iz].trms);
    
    /* 

    check for patches for this layer

    */
    for(i=0;i < nlat;i++){	/* i lat loop */
      for(j=0;j < nlon;j++){	/* j lon loop */
	
	
	os = NXY(i,j);
	if((layer[iz].patch_code[os]==0) && (temp[os] > tcut)){ /* found anomaly */
	  
	  /* should we make a new patch? */
	  /* check if neighboring patches are assigned already */
	  pfound = nbr_code(i,j,nlon,nlat,layer,iz);
	  //fprintf(stderr,"%i %i %i\n",i,j,pfound);
	  pfound = 0;
	  if(!pfound){
	    /* new patch */
	    layer[iz].npatch += 1;
	    if(layer[iz].npatch >= NMAX_PLUME){
	      fprintf(stderr,"%s: plume count out of bounds, %i\n",
		      argv[0],NMAX_PLUME);
	      exit(-1);
	    }
	    pfound = layer[iz].npatch;
	  }
	  layer[iz].patch_code[os] = pfound;
	  /* 
	     find extent of patch region 
	  */
	  midla=minla=maxla=i;
	  minlo=midlo=maxlo=j; /* start at this location */
	  //
	  for(iloop=0;iloop < 3;iloop++){
	    for(istart=minla;istart <= maxla;istart++){
	      ii = istart;
	      /* check toward the east */
	      
	      jj=j+1;if(jj==nlon)jj=0;
	      while((jj < nlon)&&((layer[iz].patch_code[NXY(ii,jj)] == 0) || 
				  (layer[iz].patch_code[NXY(ii,jj)]==layer[iz].npatch))&&
		    (temp[NXY(ii,jj)] > tcut)){
		if(layer[iz].patch_code[NXY(ii,jj)] == 0)
		  layer[iz].patch_code[NXY(ii,jj)] = layer[iz].npatch;
		if(jj>maxlo)maxlo = jj;
		if(jj<minlo)minlo = jj;
		if((maxlo == nlon-1)&&(minlo==0))
		  break;
		jj++;
		if(jj==nlon)jj=0;
	      }
	      /* check toward the west */
	      jj=minlo-1;
	      while((jj >= 0)&&((layer[iz].patch_code[NXY(ii,jj)] == 0)||
				(layer[iz].patch_code[NXY(ii,jj)] == layer[iz].npatch))&&
		    (temp[NXY(ii,jj)] > tcut)){
		if(layer[iz].patch_code[NXY(ii,jj)]==0)
		  layer[iz].patch_code[NXY(ii,jj)] = layer[iz].npatch;
		if(jj>maxlo)maxlo = jj;
		if(jj<minlo)minlo = jj;
		if((maxlo == nlon-1)&&(minlo==0))
		  break;
		jj--;
	      }
	      if(jj<0)jj=0;
	      if(jj>=nlon)jj=nlon-1;
	      
	      //fprintf(stderr,"a: %i %i %i %i\n",jj,nlon,minlo,maxlo);
	      
	      if(maxlo-minlo>.9*nlon){
		jstart = maxlo;jstop = minlo;
	      }else{
		jstart = minlo;jstop = maxlo;
	      }
	      /* done for this row */
	      for(jj=jstart;jj <= jstop;jj++){
		if(jj>=nlon)jj=0;
		ii= istart ;		/* look up */
		while((ii < nlat)&&((layer[iz].patch_code[NXY(ii,jj)] == 0)||
				    (layer[iz].patch_code[NXY(ii,jj)] == layer[iz].npatch))&&
		      (temp[NXY(ii,jj)] > tcut)){
		  if(layer[iz].patch_code[NXY(ii,jj)] == 0)
		    layer[iz].patch_code[NXY(ii,jj)] = layer[iz].npatch;
		  if(ii>maxla)maxla = ii;
		  if(ii<minla)minla = ii;
		  ii++;
		}
		ii=midla;		/* look down */
		while((ii >=0 )&&((layer[iz].patch_code[NXY(ii,jj)] == 0)||
				  (layer[iz].patch_code[NXY(ii,jj)] == layer[iz].npatch))&&
		      (temp[NXY(ii,jj)] > tcut)){
		  if(layer[iz].patch_code[NXY(ii,jj)]==0)
		    layer[iz].patch_code[NXY(ii,jj)] = layer[iz].npatch;
		  if(ii>maxla)maxla = ii;
		  if(ii<minla)minla = ii;
		  ii--;
		}
		
	      }
	      //fprintf(stderr,"b: %i %i %i %i\n",ii,nlat,minla,maxla);
	    } /* end row loop */
	      //fprintf(stderr,"c\n");
	    for(ii=minla;ii<=maxla;ii++){
	      if(maxlo-minlo>.9*nlon){
		jstart = maxlo;jstop = minlo;
	      }else{
		jstart = minlo;jstop = maxlo;
	      }
	      /* done for this row */
	      for(jj=jstart;jj <= jstop;jj++){
		if(jj>=nlon)jj=0;
		if((layer[iz].patch_code[NXY(ii,jj)] == 0) && (temp[NXY(ii,jj)] > tcut))
		  layer[iz].patch_code[NXY(ii,jj)] = layer[iz].npatch;
	      }
	    }
	  } /* end iloop */
	} /* end new plume loop */
      }	/* end lon loop */
    } /* end lat loop */
    /* 

    merge anomaly patches
    
    */	  
    for(iloop=0;iloop < 5;iloop++){
      
      for(i=0;i < nlat;i++){	/* i lat loop */
	for(j=0;j < nlon;j++){	/* j lon loop */
	  os = NXY(i,j);
	  if(temp[os] > tcut){
	    /* over limit */
	    if(!layer[iz].patch_code[os]){
	      fprintf(stderr,"%s: error: missed location %i %i\n",argv[0],i,j);
	      exit(-1);
	    }
	    
	    pfound = nbr_code(i,j,nlon,nlat,layer,iz);
	    if(pfound){		/* neighbor has been assigned */
	      if(pfound != layer[iz].patch_code[os]){
		/* need to merge */
		//fprintf(stderr,"%s: loop %i: need to merge patch %i and %i at %i,%i\n",argv[0],iloop,pfound,layer[iz].patch_code[os],i,j);
		if(pfound < layer[iz].patch_code[os]){
		  plow = pfound;phigh=layer[iz].patch_code[os];
		}else{
		  phigh= pfound;plow=layer[iz].patch_code[os];
		}
		for(ii=0;ii<nlat;ii++)
		  for(jj=0;jj<nlon;jj++)
		    if(layer[iz].patch_code[NXY(ii,jj)]==phigh)
		      layer[iz].patch_code[NXY(ii,jj)] = plow;
		    else if(layer[iz].patch_code[NXY(ii,jj)]>phigh)
		      layer[iz].patch_code[NXY(ii,jj)] = layer[iz].patch_code[NXY(ii,jj)] - 1;
		layer[iz].npatch --;
	      }
	    }
	  }
	}
      }
    }

    /* 

    count plume samples, their area, and compute centroids in i,j space

    */
    layer[iz].parea     = (GGRD_CPREC *)calloc((size_t)layer[iz].npatch,sizeof(GGRD_CPREC));
    layer[iz].ppixcount = (unsigned int *)calloc(layer[iz].npatch,sizeof(unsigned int));
    layer[iz].pcentroid = (struct pxl *)calloc(layer[iz].npatch,sizeof(struct pxl));
    xp = (GGRD_CPREC *)realloc(xp,3*layer[iz].npatch*sizeof(GGRD_CPREC));
    if(!xp || !layer[iz].parea || !layer[iz].ppixcount || !layer[iz].pcentroid ){
      fprintf(stderr,"%s: mem error\n",argv[0]);exit(-1);}
    for(i=0;i<layer[iz].npatch;i++){ /* for centroids */
      for(j=0;j<3;j++)
	xp[i*3+j]=0.0;
    }
    r = nd_r_of_zkm(layer[iz].z);
    r2 = r * r;
    for(i=os=0;i<nlat;i++){
      tmp = w[i] * r2;
      for(j=0;j<nlon;j++,os++){
	if(layer[iz].patch_code[os]){
	  //fprintf(stdout,"%6g %.2f %4i %g\n",lon[j],lat[i],layer[iz].patch_code[os],w[i]);
	  apatch = layer[iz].patch_code[os] - 1; /* patch code */
	  layer[iz].ppixcount[apatch] += 1;
	  layer[iz].parea[apatch] += tmp; /* area per patch */
	  layer[iz].pa += tmp;		/* total */
	  lonlat2xyz(lon[j],lat[i],xloc);
	  for(k=0;k<3;k++)	/* add for centroid */
	    xp[apatch*3+k] += xloc[k];
	}
      }
    }
    /* total patch area normalized by surface (not shell) area */
    layer[iz].pa /= ws;

    for(i=0;i < layer[iz].npatch;i++){
      layer[iz].parea[i] /= ws;
      /* compute lonlat */
      if(layer[iz].ppixcount[i]){
	for(k=0;k < 3;k++)
	  xp[i*3+k] /= (double)layer[iz].ppixcount[i];
	xyz2lonlat((xp+i*3),&tmplon,&tmplat);
	/* 
	   convert centroid into i,j space 
	*/
	layer[iz].pcentroid[i].c[0]   = tmplon/dlon; /* lon of centroid */
	layer[iz].pcentroid[i].c[1] = (tmplat-lat[0])/dlat;
	//fprintf(stderr,"%3i %11g %11g %11g %11g\n",i+1,tmplon,tmplat,lon[layer[iz].pcentroid[i].c[0]],lat[layer[iz].pcentroid[i].c[1]]);
      }
    }
    /*  */
    fprintf(stderr,"z: %6g pfrac: %6.2f%% nrp: %5i detected tmean: %3.f tmax: %.3f dt: %.3f trms: %.3f tcut: %.3f\n",
	    layer[iz].z,layer[iz].pa*100,
	    layer[iz].npatch,layer[iz].tmean,
	    layer[iz].tmax,layer[iz].tmax-layer[iz].tmean,
	    layer[iz].trms,tcut);
    if(0){
      if(iz == 0)
	out = fopen("patch.dat","w");
      /* for debugging, print all patches */
      print_patch((layer+iz),out,lon,lat,nlon,nlat);
      if(iz==nz-1)
	fclose(out);
    }

  } /* end depth loop */
 
  /* depth range */
  zrange = layer[0].z - layer[nz-1].z;
  
  /* weight for each layer in volume summation */
  dzw = 1.0 / (float)nz;

  fprintf(stderr,"%s: patch detection done, now assigning plumes within %g km shell\n",
	  argv[0],zrange);
  /* 

  look for continuous plumes

  */
  if(upward){
    /* 
       
    search upward, old method

    */
    plume = NULL;
    nplume = 0;
    for(iz=0;iz < nz;iz++){		
      /* 
	 depth loop 
      */
      for(k=0;k < layer[iz].npatch;k++){	/* patches per depth layer loop */
	apatch = k + 1;
	nfound = 0;
	found = old_found = -1;
	if(iz > 0){
	  /* 
	     check if this patches' location is found in any plume of
	     the layer below
	  */
	  
	  old_found = -1;
	  loop = 1;
	  while(loop){
	    loop = 0;
	    for(nfound = i = 0;i < nplume;i++){
	      /* check if patch apatch in layer iz is found in any
		 plumes layer below */
	      found = check_for_patch(i,plume,(iz-1),iz,apatch,layer,
				      nlon,nlat,rmax,lon,lat,layer[iz].z);
	      
	      if(found != -1){
		nfound++;
		if(!allow_downward_branch){ /* exit if we're not allowing
					       downward branching */
		  old_found = found;
		  break;
		}
		if((old_found != -1)&&(old_found != found)){
		  merge_plumes(&plume,&nplume,old_found,found,layer,nlon,nlat);
		  loop = 1;
		  break;
		}else{
		  old_found = found;
		}
	      }
	    }
	  } /* loop */
	}
	if(nfound == 0){
	  /* 
	     make new plume
	  */
	  make_new_plume(&plume,&nplume);
	  /* add this patch to the new plume */
	  //fprintf(stderr,"adding %i to new p %i at level %i\n",apatch,nplume,iz);
	  add_patch_to_plume(iz,(nplume-1),apatch,layer,plume,nlon,nlat);	/*  */
	  
	}else{
	  /* use old one as stored in found */
	  //fprintf(stderr,"adding %i to old p %i at level %i\n",apatch,old_found,iz);
	  add_patch_to_plume(iz,old_found,apatch,layer,plume,nlon,nlat); /*  */
	}
      } /* end k (number of patch per layer loop)*/
    } /* end iz loop */
  }else{
    /* 
       downward branch, allow multiple assignments of patches to
       plumes
    */
    plume = NULL;
    nplume = 0;
    for(iz=nz-1;iz >= 0;iz--){		
      /* 
	 depth loop 
      */
      for(k=0;k < layer[iz].npatch;k++){	/* patches per depth layer loop */
	apatch = k + 1;

	nfound = 0;
	if(iz < nz-1){
	  /* 
	     check if this patches' location is found in any plume of
	     the layer above
	  */
	  for(i = 0;i < nplume;i++){
	    found =  check_for_patch(i,plume,(iz+1),iz,apatch,layer,
				     nlon,nlat,rmax,lon,lat,layer[iz].z);
	    if(found != -1){
	      add_patch_to_plume(iz,found,apatch,layer,plume,nlon,nlat); /*  */
	      nfound++;
	    }
	  }
	}
	if(!nfound){
	  /* 
	     make new plume
	  */
	  make_new_plume(&plume,&nplume);
	  /* add this patch to the new plume */
	  add_patch_to_plume(iz,(nplume-1),apatch,layer,plume,nlon,nlat);	/*  */
	}
	//fprintf(stderr,"%i %g %i %i\n",iz,layer[iz].z,apatch,nplume);
      } /* end k (number of patch per layer loop)*/
    } /* end iz loop */
  }

  /* 
     compute mean properties and sort plumes
  */
  compute_mean_plume_properties(plume,nplume,w,ws,layer,dzw);
  
  /* 
     output all non-zero length plumes to a list  
  */
  out2 = fopen("plumelist.dat","w");
  fprintf(out2,"# plume_nr zmin zmax length total_area mean_area total_volume\n");
  nzero = 0;
  for(i=0;i < nplume;i++)
    if(plume[i].zmax - plume[i].zmin > 0){
      nzero++;
      fprintf(out2,"%5i\t %11g %11g %11g\t  %8.3e  %8.3e %8.3e\n",
	      i+1,
	      plume[i].zmin,plume[i].zmax,(plume[i].zmax-plume[i].zmin),
	      plume[i].total_area,plume[i].mean_area, plume[i].volume);
    }
  fclose(out2);
  fprintf(stderr,"%s: written all %i non-zero length plumes to plumelist.dat\n",
	  argv[0],nzero);
  /* 

  print the biggest plumes to separate files
  
  */
  fprintf(stderr,"%s: output if length > %g and zmin < %g\n",
	  argv[0],min_length,max_depth);
  nplume_select=0;
  tvol = marea = 0.0;
  passigned=(PCTYPE *)calloc((size_t)(nn*nz),sizeof(PCTYPE));
  if(!passigned){fprintf(stderr,"%s: mem error\n",argv[0]);exit(-1);}

  /* remove old plume dat files */
  system("rm plume.???.dat");
  for(i=0;i<nplume;i++){
    if(((plume[i].zmax - plume[i].zmin) >= min_length)&&(plume[i].zmin <= max_depth)){
      /* 
	 write plume pixels to file 
      */
      sprintf(filename,"plume.%03i.dat",nplume_select+1);
      out = fopen(filename,"w");
      fprintf(stderr,"%s: printing plume %3i to %s, depth extent: %10.1f - %10.1fkm, vol: %g\n",
	      argv[0],nplume_select+1,filename,
	      plume[i].zmin,plume[i].zmax,
	      plume[i].volume);
      print_plume((plume+i),out,lon,lat,layer);
      fclose(out);
      /* sum up */
      marea += plume[i].mean_area;
      tvol += plume[i].volume;
      nplume_select++;
      /* 
	 check for pixel assignment
      */
      add_plume_to_px_assign((plume+i),nlon,nlat,layer,passigned);
    }
  }
  marea /= (double)nplume_select;
  /* compute the fractional volume of the assigned pixels */
  rvol = compute_total_pix_assigned_vol(nlon,nlat,nz,layer,passigned,w,ws,dzw);
  /* 

  main output

  */
  fprintf(stderr,"%s: nr long plume: %i out of: %i mean area: %g%% of Earth's surface, volume: %g%%, corr. col: %g%% \n",
	  argv[0],nplume_select,nplume,marea*100., tvol * 100,
	  rvol * 100);

  /* output of 
     
     number of selected plumes

     mean plume area, fractional of selected plumes
     
     fractional volume of selected plumes, possibly double counting 

     total number of plumes, including short ones

     real plume volume, fractional, of selected plumes

  */
  out = fopen("plume.dat","w");
  fprintf(out,"%i %g %g %i %g\n",nplume_select,marea,tvol,nplume,rvol);
  fclose(out);
  exit(0);
}

/* merge two plumes */

void merge_plumes(struct plm **plume, int *nplume,
		  int aplume, int bplume,
		  struct lyr *layer, int nlon, int nlat)
{
  struct plm *use_plume;
  int i,npn,j,k;
  fprintf(stderr,"merge_plumes: merging plumes %i and %i, out of total of %i\n",
	  aplume+1,bplume+1,*nplume);
  
  /* logic check */
  if((aplume > *nplume-1)||(bplume > *nplume-1)||(aplume < 0)||(bplume <0)){
    fprintf(stderr,"merge_plumes: a: %i b: %i out of bounds, only %i total plumes\n",
	    aplume,bplume,*nplume);
    exit(-1);
  }
  /* make room for a new plume */
  make_new_plume(plume,nplume);

  /* add first plume */
  use_plume = (*plume+aplume);
  for(i=0;i < use_plume->nz;i++) /* loop through all levels in this plume */
    for(j=0;j < use_plume->pnr[i];j++){ /* loop through all patches
					   selected for this layer */
      //fprintf(stderr,"o1: %i %i %i\n",aplume,i,use_plume->pname[i][j]);
      add_patch_to_plume(use_plume->iz[i],(*nplume-1),
			 use_plume->pname[i][j],
			 layer,*plume,nlon,nlat);	/*  */

    }
  /* second plume */
  use_plume = (*plume+bplume);
  for(i=0;i < use_plume->nz;i++) /* loop through all levels in this plume */
    for(j=0;j < use_plume->pnr[i];j++){ /* loop through all patches
					   selected for this layer */
      //fprintf(stderr,"o2: %i %i %i\n",bplume,i,use_plume->pname[i][j]);
      add_patch_to_plume(use_plume->iz[i],(*nplume-1),
			 use_plume->pname[i][j],
			 layer,*plume,nlon,nlat);	/*  */
    }
  npn = 0;
  for(i=0;i< *nplume;i++)
    if((i != aplume) && (i != bplume)){
      if(i > npn){
	/* make room and copy down */
	use_plume = (*plume+i);
	/* clear old */
	clear_plume((*plume + npn));
	/* overwrite with this one */
	for(j=0;j < use_plume->nz;j++) /* loop through all levels in this plume */
	  for(k=0;k < use_plume->pnr[j];k++) /* loop through all patches
					       selected for this layer */
	    add_patch_to_plume(use_plume->iz[j],npn,
			       use_plume->pname[j][k],
			       layer,*plume,nlon,nlat);	
      }
      npn++;
    }
  /* free the remainders */
  for(i=npn;i < *nplume;i++)
    clear_plume((*plume+i));
  *nplume = npn;
  fprintf(stderr,"merge_plumes: now %i plumes\n",*nplume);
  if(0){
    /* print out combined */
    use_plume = (*plume+(*nplume - 1));
    for(i=0;i < use_plume->nz;i++) /* loop through all levels in this plume */
      for(j=0;j < use_plume->pnr[i];j++){ /* loop through all patches
					     selected for this layer */
	fprintf(stderr,"n: %i %i\n",i,use_plume->pname[i][j]);
      }
  }
}
/* 

determine if plume iplume has any patches at a location roughly rmax
distance from the patch code apatch in layer ilayer and the plume at
layer ilayerp

will return the number of the plume as passsed (iplume) if plume
found, -1 else

*/
int check_for_patch(int iplume, struct plm *plume,
		    int ilayerp,int ilayer, 
		    int apatch, struct lyr *layer, 
		    int nlon, int nlat,
		    GGRD_CPREC rmax,GGRD_CPREC *lon, 
		    GGRD_CPREC *lat,GGRD_CPREC z)
  {
  int i,j,k,have_layer;

  /* check if this plume has any entries in the layer  ilayer */
  have_layer = -1;
  for(i=0;i < plume[iplume].nz;i++)
    if(plume[iplume].iz[i] == ilayerp){
      have_layer = i;break;
    }
  if(have_layer != -1){
    /* 
       the layer is in plume, check if the patch location is
       represented
    */
    for(i=0;i<nlat;i++){	/* lat loop */
      for(j=0;j<nlon;j++){	/* lon loop */
	if(layer[ilayer].patch_code[NXY(i,j)] == apatch){
	  for(k=0;k<plume[iplume].npixel[have_layer];k++){
	    if(int_dist(plume[iplume].pixels[have_layer][k].c[1],
			plume[iplume].pixels[have_layer][k].c[0],
			i,j,lon,lat,z) < rmax)
	      return iplume;
	  }
	}	  
      }
    }
  }
  return -1;			/* not found */
}
/* 

add all pixels within the patch apatch to this plume


ilayer = current layer
iplume = current plume
apatch = patch code
layer = layer structure
plume = plume structure 


*/
void add_patch_to_plume(int ilayer, int iplume,
			int apatch, struct lyr *layer, 
			struct plm *plume, int nlon,
			int nlat)
{
  int i,j,luse,puse;
  static int print = 0;
  luse = -1;
  for(i=0;i < plume[iplume].nz;i++){
    if(plume[iplume].iz[i] == ilayer){
      luse = i;break;
    }
  }
  if(luse == -1){
    /* 
       add new layer code to list 
    */
    plume[iplume].iz = 
      (int *)realloc(plume[iplume].iz,(plume[iplume].nz+1)*sizeof(int));
    plume[iplume].npixel = 
      (int *)realloc(plume[iplume].npixel,(plume[iplume].nz+1)*sizeof(int));
    /*  */
    plume[iplume].pnr = 
      (int *)realloc(plume[iplume].pnr,(plume[iplume].nz+1)*sizeof(int));
    plume[iplume].pname = 
      (int **)realloc(plume[iplume].pname,(plume[iplume].nz+1)*sizeof(int *));
    /*  */
    luse = plume[iplume].nz;
    
    plume[iplume].pnr[luse] = 0; /* number of different patches per layer */
    plume[iplume].npixel[luse] = 0;
    plume[iplume].iz[luse] = ilayer;
    /* add new pixel point */
    plume[iplume].pixels = 
      (struct pxl **)realloc(plume[iplume].pixels,sizeof(struct pxl *)*(plume[iplume].nz+1));
    if(!plume[iplume].iz || !plume[iplume].npixel || !plume[iplume].pixels || !plume[iplume].pnr){
      fprintf(stderr,"mem error \n");exit(-1);}
    plume[iplume].pixels[luse] = NULL;
    /*  */
    plume[iplume].nz++;  

  }
  /* check if we have this patch number in the list already */
  for(i=0;i < plume[iplume].pnr[luse];i++){
    if(plume[iplume].pname[luse][i] == apatch){
      puse = i;break;
    }
  }
  if(i == plume[iplume].pnr[luse]){ 
    /* 
       didn't find this patch code, add to list 
    */
    puse = plume[iplume].pnr[luse];
    if(puse == 0)
      plume[iplume].pname[luse] = (int *)malloc(sizeof(int));
    else
      plume[iplume].pname[luse] = (int *)realloc(plume[iplume].pname[luse],
						 sizeof(int)*(plume[iplume].pnr[luse]+1));
    plume[iplume].pname[luse][puse] = apatch;
    plume[iplume].pnr[luse] += 1;
  }
  if(print)
    fprintf(stderr,"add pat %3i layer %3i to plume %3i, local layer %3i(%3i), patches %3i(%3i) npix_in: %5i ",
	    apatch,ilayer+1,iplume+1,
	    luse+1,plume[iplume].nz,
	    puse+1,plume[iplume].pnr[luse],
	    plume[iplume].npixel[luse]);
  for(i=0;i<nlat;i++){	/* lat loop */
    for(j=0;j<nlon;j++){	/* lon loop */
      if(layer[ilayer].patch_code[NXY(i,j)] == apatch){
	/* add pixels to list */
	plume[iplume].pixels[luse] = (struct pxl *)
	  realloc(plume[iplume].pixels[luse],sizeof(struct pxl)*(plume[iplume].npixel[luse]+1));
	plume[iplume].pixels[luse][plume[iplume].npixel[luse]].c[0] = j; /* lon */
	plume[iplume].pixels[luse][plume[iplume].npixel[luse]].c[1] = i; /* lat */
	plume[iplume].npixel[luse]++;
      }
    }
  }
  if(print)
    fprintf(stderr," nout: %5i\n",plume[iplume].npixel[luse]);
}

void print_plume(struct plm *plume,FILE *out,GGRD_CPREC *lon, GGRD_CPREC *lat,struct lyr *layer)
{
  int i,j;
  for(i=0;i<plume->nz;i++){
    for(j=0;j<plume->npixel[i];j++)
      fprintf(out,"%11g %11g %11g\n",
	      lon[plume->pixels[i][j].c[0]],lat[plume->pixels[i][j].c[1]],
	      layer[plume->iz[i]].z);
    fprintf(out,"\n");
  }
  
}
/* print patches of this layer */
void print_patch(struct lyr *layer,FILE *out,GGRD_CPREC *lon, GGRD_CPREC *lat,int nlon, int nlat)
{
  int i,j,os;
  for(i=os=0;i<nlat;i++){
    for(j=0;j<nlon;j++,os++){
      if(layer->patch_code[os]){
	fprintf(out,"%11g %11g %11g\n",
		lon[j],lat[i],layer->z);
      }
    }
  }
  fprintf(out,"\n");
}

void add_plume_to_px_assign(struct plm *plume,int nlon, int nlat, struct lyr *layer,
			    PCTYPE *passigned)
{
  int i,j,iz,ioff;
  static int called = 0,nn;
  if(!called){
    nn = nlon*nlat;
    called = 1;
  }
  for(i=0;i<plume->nz;i++){
    iz = plume->iz[i];
    for(j=0;j<plume->npixel[i];j++){
      ioff = iz * nn + 
	plume->pixels[i][j].c[0] * nlat + 
	plume->pixels[i][j].c[1];
      if(!passigned[ioff])
	passigned[ioff] = 1;
    }
  }
}


GGRD_CPREC compute_total_pix_assigned_vol(int nlon,int nlat,int nz,struct lyr *layer,
					  PCTYPE *passigned,GGRD_CPREC *w, GGRD_CPREC ws,
					  GGRD_CPREC dz)
{
  GGRD_CPREC rvol,r,r2,rarea;
  int i,j,k,nn,off;
  nn = nlat * nlon;
  rvol = 0;
  for(i=0;i<nz;i++){
    r =  nd_r_of_zkm(layer[i].z); /* radius */
    r2 = r*r;			/* r^2 */
    rarea = 0.0;
    for(j=0;j<nlon;j++)
      for(k=0;k<nlat;k++){
	off = i * nn + j*nlat + k;
	if(passigned[off]){
	  rarea += w[k];
	}
      }
    rarea = rarea *r2/ws;
    rvol += rarea*dz;
  }

  return rvol;


}
/* 

compute the mean conduit properties for all plumes 

*/
void compute_mean_plume_properties(struct plm *plume, int nplume,
				   GGRD_CPREC *w, GGRD_CPREC ws,
				   struct lyr *layer,GGRD_CPREC dz)
{
  int i,j,k;
  GGRD_CPREC r,z;
  for(i=0;i < nplume;i++){	/* loop through plumes */
    /* make room for area computation */
    plume[i].area = (GGRD_CPREC *)
      realloc(plume[i].area,sizeof(GGRD_CPREC)*plume[i].nz);
    if(!plume[i].area){fprintf(stderr,"mem error\n");exit(-1);}
    plume[i].zmin = 1e20;plume[i].zmax = -1e20;

    /* this shouldn't be necessary, but hey */
    plume[i].mean_npixel = 0;
    plume[i].volume = 0.0;
    plume[i].total_area = 0.0;
    
    for(j=0;j < plume[i].nz;j++){	/* loop through layers */

      z = layer[plume[i].iz[j]].z; /* this layer depth */
      /* add number of pixels */
      plume[i].mean_npixel += plume[i].npixel[j];
      /* sum up all pixel contributions for this level */
      plume[i].area[j] = 0.0;
      for(k=0;k < plume[i].npixel[j];k++){ /* add up all pixels for area */
	//fprintf(stderr,"%i %i %i %i %g\n",i,j,k,plume[i].pixels[j][k].c[1],w[plume[i].pixels[j][k].c[1]]);
	plume[i].area[j] += w[plume[i].pixels[j][k].c[1]];
      }
      /* 
	 area of plume entries for this layer 
      */
      r = nd_r_of_zkm(z);
      
      plume[i].area[j] = plume[i].area[j]/ws * r*r; /* scale with r^2 */

      /* total plume area */
      plume[i].total_area += plume[i].area[j];
      /* total volume */
      plume[i].volume += plume[i].area[j] * dz;
            
      //fprintf(stderr,"%i %i %g %g %g %i\n",i,j,z,r,plume[i].area[j]*100,plume[i].npixel[j]);
      
      if(z > plume[i].zmax)
	plume[i].zmax = z;
      if(z < plume[i].zmin)
	plume[i].zmin = z;
    }
    /* mean number of pixels */
    plume[i].mean_npixel = 
      (int)((float)plume[i].mean_npixel/(float)plume[i].nz+.5);
    /* mean plume area */
    plume[i].mean_area = 
      plume[i].total_area / (float)plume[i].nz;
  }

  /* sort plumes by their maximum number of entries */
  qsort(plume,nplume,sizeof(struct plm),plm_compare);

}
int plm_compare(struct plm *a, struct plm *b)
{
  if(a->nz > b->nz)
    return -1;
  else if(a->nz == b->nz)
    return 0;
  else 
    return 1;

}
void make_new_plume(struct plm **plume,int *nplume)
{
  *plume = (struct plm *)realloc(*plume, (*nplume+1)*sizeof(struct plm));

  if(!(*plume)){
    fprintf(stderr,"mem error plume %i\n",*nplume+1);
    exit(-1);
  }
  init_plume((*plume+ *nplume));
  *nplume += 1;
}
/* clear a single plume structure */
void clear_plume(struct plm *plume)
{
  int i;
  free(plume->iz);
  for(i=0;i< plume->nz;i++){
    free(plume->pname[i]);
    free(plume->pixels[i]);
  }
  free(plume->pname);
  free(plume->pixels);
  free(plume->npixel);
  free(plume->pnr);
  free(plume->area);
  init_plume(plume);
}
/* get rid of single plume  */
void free_plume(struct plm **plume)
{
  clear_plume(*plume);
  free(*plume);
  
}

void init_plume(struct plm *plume)
{
  (plume)->nz = 0;
  (plume)->iz = NULL;
  (plume)->pixels = NULL;
  (plume)->npixel = NULL;
  (plume)->pname = NULL;
  (plume)->pnr = NULL;
  (plume)->area = NULL;
  (plume)->mean_area = 0.0;
  (plume)->mean_npixel = 0;
}

/* give z in km, lon and lat arrays in degrees

output distance will be in km */
GGRD_CPREC int_dist(int i1,int j1,int i2,int j2,
	       GGRD_CPREC *lon,GGRD_CPREC *lat,
	       GGRD_CPREC z)
{
  GGRD_CPREC lon1,lat1,lon2,lat2,r,tmp1,tmp2,tmp3;
  static GGRD_CPREC pif = GGRD_PI/180;
  r = nd_r_of_zkm(z) * R_EARTH;		
  
  lon1=lon[j1]*pif;
  lon2=lon[j2]*pif;
  lat1=lat[i1]*pif;
  lat2=lat[i2]*pif;

  tmp1=sin((lat1-lat2)/2.0);
  tmp1=tmp1*tmp1;
  tmp2=sin((lon1-lon2)/2.0);
  tmp2=tmp2*tmp2;
  tmp2*=cos(lat1);
  tmp2*=cos(lat2);
  tmp3=sqrt(tmp1+tmp2);
  return r * 2.0*asin(tmp3);
}

/* convert lon lat in degrees into cartesian */
void lonlat2xyz(GGRD_CPREC lon, GGRD_CPREC lat,
		GGRD_CPREC *xc)
  {
  GGRD_CPREC lambda,phi,tmp;
  static double pif = GGRD_PI/180;
  lambda=lat*pif;
  phi=lon*pif;
  
  tmp=cos(lambda);
    
  xc[0]=tmp * cos(phi);
  xc[1]=tmp * sin(phi);
  xc[2]=sin(lambda);
}
/*  */
void xyz2lonlat(GGRD_CPREC *xc,GGRD_CPREC *lon, GGRD_CPREC *lat)
{
  GGRD_CPREC tmp1,phi,theta;
  static double f = 180/GGRD_PI,twopi=GGRD_PI*2.;

  tmp1 = xc[0]*xc[0] + xc[1]*xc[1];
  /*  */
  theta=atan2(sqrt(tmp1),xc[2]);
  phi=atan2(xc[1],xc[0]);
  if(phi < 0)
    phi += twopi;
  *lon = phi * f;
  *lat = 90. - theta * f;
}
GGRD_CPREC nd_r_of_zkm(GGRD_CPREC z)
{
  return 1.0-fabs(z)/R_EARTH;
}
/* 
   check for neighboring patch anomalies 
*/
PCTYPE nbr_code(int i, int j, int nlon, int nlat, struct lyr *layer, int iz)
{
  int os,left,right,top,down;
  GGRD_CPREC r;
  PCTYPE p;
  if(j==0)left = nlon-1;else left = j-1;
  if(j==nlon-1)right =0;else right= j+1;
  //  fprintf(stderr,"%i %i %i %i\n",layer[iz].patch_code[NXY(i,left)],layer[iz].patch_code[NXY(i,right)],(i>1)?(layer[iz].patch_code[NXY(i-1,j)]):(0),(i<nlat-1)?(layer[iz].patch_code[NXY(i+1,j)]):(0));
  if((p = layer[iz].patch_code[NXY(i,left)]))
    return p;
  if((p = layer[iz].patch_code[NXY(i,right)]))
    return p;
  if(i > 0){
    if((p = layer[iz].patch_code[NXY(i-1,j)]))
      return p;
    if(0){
      if((p = layer[iz].patch_code[NXY(i-1,left)]))
	return p;
      if((p = layer[iz].patch_code[NXY(i-1,right)]))
	return p;
    }
  }
  if(i<nlat-1){
    if((p = layer[iz].patch_code[NXY(i+1,j)]))
      return p;
    if(0){
      if((p = layer[iz].patch_code[NXY(i+1,left)]))
	return p;
      if((p = layer[iz].patch_code[NXY(i+1,right)]))
	return p;
    }
  }

  return 0;
  
}
