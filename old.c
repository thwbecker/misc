/* 

using GMT grd files in 3-D for tomographic models, detect
plumes 

$Id$

*/
#include "ggrd_grdtrack_util.h"
#include <limits.h>

#define NXY(i,j) (i*nlon + j)

/* use unsigned char */
#define PCTYPE unsigned char
#define NMAX_PLUME UCHAR_MAX



/* pixcel, with two integer coord */
struct pxl{
  unsigned int c[2];		/* j (0...nlat-1), i (0...nlon-1)  */
};
/* 
   plume structure 
*/
struct plm{
  int nz,*iz,			/* number of depth levels, and depth level code */
    *npixel;			/* number of pixels per layer */
  struct pxl **pixels;		/* pixels[nz] */
  int *pnr,**pname;		/* patch numbers and names */
  int mean_npixel;
  GGRD_CPREC *area,mean_area,zmin,zmax;
};
/* layer structure  */
struct lyr{
  GGRD_CPREC z,pa;			/* depth, total patch area*/
  GGRD_CPREC tmean,tmax;	/* extrema */
  PCTYPE *patch_code;	/* all codes, which ones assigned */
  int npatch;
  /* per patch */
  unsigned int *ppixcount;
  struct pxl *pcentroid;
  GGRD_CPREC *parea;
};

int check_for_patch(int , struct plm *,
		    int ,int , struct lyr *, 
		    int , int ,GGRD_CPREC ,GGRD_CPREC *, 
		    GGRD_CPREC *,GGRD_CPREC );
void lonlat2xyz(GGRD_CPREC , GGRD_CPREC ,GGRD_CPREC *);
void xyz2lonlat(GGRD_CPREC *,GGRD_CPREC *, GGRD_CPREC *);
void make_new_plume(struct plm **,int *);
void add_patch_to_plume(int , int , int,struct lyr *, struct plm *, int ,
			int );
GGRD_CPREC int_dist(int ,int ,int ,int ,GGRD_CPREC *,GGRD_CPREC *,
	       GGRD_CPREC );
void compute_mean_plume_properties(struct plm *, int,GGRD_CPREC *, GGRD_CPREC ,struct lyr *);
void print_plume(struct plm * ,FILE *,GGRD_CPREC *, GGRD_CPREC *,struct lyr *);
int plm_compare(struct plm *, struct plm *);


int main(int argc, char **argv)
{
  GGRD_CPREC *lon,*lat,ws,*w,frac,*temp,dlon,dlat,dt,tcut,
    dz,*xp=NULL,xloc[3],tmplon,tmplat,rmax,tvol,marea,dfac;
  struct ggrd_gt t3d[1];
  double value;
  static double pif = GGRD_PI/180, scale = 1.0;	/* temp/vel anom scale */
  int i,j,k,nlon,nlat,os,iz,nn,ii,jj,minlo,maxlo,minla,maxla,found,
    midlo,midla,istart,iloop,jstart,jstop,jloop,nz,apatch,nplume;
  struct lyr *layer;
  struct plm *plume = NULL;
  char filename[300];
  FILE *out;
  static GGRD_CPREC zbot = 2800, ztop = 300;
  if(argc < 3){
    fprintf(stderr,"usage:\n\n%s prefix depth_file\n",argv[0]);
    
    exit(-1);
  }
  /* 
     plume detection settings 
  */
  frac = .5;			/* cutoff level  */

  dz = 50;			/* depth spacing */

  dfac = 3;			/* how many depth spacings to allow
				   vertically? */



  /* parameterization settings */
  dlon = dlat = 1;		/* spacing of interpolation grid */


  if(argc > 3)
    sscanf(argv[3],"%lf",&frac);

  fprintf(stderr,"%s: prefix: %s dfile: %s, using frac: %g\n",
	  argv[0],argv[1],argv[2],frac);
  
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
  if(ggrd_grdtrack_init_general(TRUE,argv[1], argv[2],"-L",t3d,FALSE,FALSE)){
    fprintf(stderr,"%s: ggrd init error\n",argv[0]);
    exit(-1);
  }
  /* 
     loop through layers bottom up
  */
  for(iz=0;iz < nz;iz++){
  //for(iz=0;iz < 1;iz++){
    if(iz==0)
      layer[iz].z = zbot;
    else
      layer[iz].z = layer[iz-1].z - dz;
    /* 
       allocate memory for this layer and set to zero
    */
    layer[iz].patch_code = (PCTYPE *)calloc(nn,sizeof(PCTYPE));
    if(!layer[iz].patch_code){fprintf(stderr,"%s: mem error\n",argv[0]);exit(-1);}
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
	temp[os] = (1 - value/100) * scale;
	/* layer mean */
	layer[iz].tmean += w[i] * temp[os];
	/* layer max */
	if(layer[iz].tmax < temp[os])
	  layer[iz].tmax = temp[os];
      }
    }
    layer[iz].tmean /= ws;
    /* cutoff temperatures */
    dt = layer[iz].tmax - layer[iz].tmean;
    tcut = frac*dt + layer[iz].tmean;	/* cutoff temperature */
    //fprintf(stderr,"%s: z: %g tmean: %g tmax: %g\n",argv[0],layer[iz].z,layer[iz].tmean,layer[iz].tmax);
    /* 

    check for plumes for this layer

    */
    for(i=0;i<nlat;i++){	/* lat loop */
      for(j=0;j<nlon;j++){	/* lon loop */


	os = NXY(i,j);
	if((layer[iz].patch_code[os]==0) && (temp[os] > tcut)){ /* found plume */
	  /* new plume */
	  layer[iz].npatch += 1;
	  if(layer[iz].npatch >= NMAX_PLUME){
	    fprintf(stderr,"%s: plume count out of bounds, %i\n",
		    argv[0],NMAX_PLUME);
	    exit(-1);
	  }
	  layer[iz].patch_code[os] = layer[iz].npatch;
	  /* 
	     find extent of plume region 
	  */
	  midla=minla=maxla=i;
	  minlo=midlo=maxlo=j; /* start at this location */
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
		jj++;if(jj==nlon)jj=0;
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
		jj--;
	      }
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
	    } /* end row loop */
	    for(ii=minla;ii<=maxla;ii++)
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
	} /* end new plume loop */
      }	/* end lon loop */
    } /* end lat loop */
    fprintf(stderr,"%s: patch detection done, now assigning plumes\n",
	    argv[0]);
    /* 

    count plume samples, their area, and compute centroids in i,j space

    */
    layer[iz].parea     = (GGRD_CPREC *)calloc(layer[iz].npatch,sizeof(GGRD_CPREC));
    layer[iz].ppixcount = (unsigned int *)calloc(layer[iz].npatch,sizeof(unsigned int));
    layer[iz].pcentroid = (struct pxl *)calloc(layer[iz].npatch,sizeof(struct pxl));
    xp = (GGRD_CPREC *)realloc(xp,3*layer[iz].npatch*sizeof(GGRD_CPREC));
    if(!xp || !layer[iz].parea || !layer[iz].ppixcount || !layer[iz].pcentroid ){
      fprintf(stderr,"%s: mem error\n",argv[0]);exit(-1);}
    for(i=0;i<layer[iz].npatch;i++){ /* for centroids */
      for(j=0;j<3;j++)
	xp[i*3+j]=0.0;
    }
    for(i=os=0;i<nlat;i++)
      for(j=0;j<nlon;j++,os++){
	//fprintf(stdout,"%6g %.2f %4i\n",lon[j],lat[i],layer[iz].patch_code[os]);
	if(layer[iz].patch_code[os]){
	  apatch = layer[iz].patch_code[os] - 1; /* plume code */
	  layer[iz].ppixcount[apatch] += 1;
	  layer[iz].parea[apatch] += w[i]; /* area per plume */
	  layer[iz].pa += w[i];		/* total */
	  lonlat2xyz(lon[j],lat[i],xloc);
	  for(k=0;k<3;k++)	/* add for centroid */
	    xp[apatch*3+k] += xloc[k];
	}
      }

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
    fprintf(stderr,"z: %6g pfrac: %.2f%% nrp: %5i detected tmean: %3.f tmax: %.3f\n",
	    layer[iz].z,layer[iz].pa*100,layer[iz].npatch,layer[iz].tmean,layer[iz].tmax);
  } /* end depth loop */

  /* 

  look for continuous plumes without allowing for branching

  */
  plume = NULL;
  nplume = 0;
  for(iz=0;iz < nz;iz++){		
    /* depth loop */
    for(k=0;k < layer[iz].npatch;k++){	/* patches per depth layer loop */
      apatch = k + 1;
      found = -1;
      if(iz > 0){
	/* check if this patches' location is found in plume check
	   layer below */
	for(i=0;i < nplume;i++){
	  found = check_for_patch(i,plume,iz, apatch,layer,
				  nlon,nlat,rmax,
				  lon,lat,layer[iz].z);
	  if(found != -1)
	    break;
	}
      }
      if(found == -1){
	/* make new plume*/
	make_new_plume(&plume,&nplume);
	add_patch_to_plume(iz,(nplume-1),apatch,layer,plume,nlon,nlat);	/*  */
      }else{
	add_patch_to_plume(iz,found,apatch,layer,plume,nlon,nlat); /*  */
      }
    }
  }
  /* 
     compute mean properties and sort plumes
  */
  compute_mean_plume_properties(plume,nplume,w,ws,layer);
  /* 
     print the biggest plumes
  */
  i=0;tvol = marea = 0.0;
  while(plume[i].nz*dz > 500){
    sprintf(filename,"plume.%03i.dat",i+1);
    out = fopen(filename,"w");
    tvol += (plume[i].zmax-plume[i].zmin)*plume[i].mean_area;
    fprintf(stderr,"%s: printing plume %3i to %s, depth extent: %10.1f - %10.1fkm, vol: %g\n",
	    argv[0],i+1,filename,plume[i].zmin,plume[i].zmax,
	    (plume[i].zmax-plume[i].zmin)*plume[i].mean_area);
    print_plume((plume+i),out,lon,lat,layer);
    fclose(out);
    marea += plume[i].mean_area;
    i++;
  }
  marea /= (double)i;
  /* 

  main output

  */
  fprintf(stderr,"%s: nr long plume: %i area: %g volume: %g \n",argv[0],i,marea, tvol);
  out = fopen("plume.dat","w");
  fprintf(out,"%i %g %g\n",i,marea,tvol);
  fclose(out);
  exit(0);
}
/* 

determine if this plume has any patches at a location roughly rmax
distance from the patch code apatch in layer ilayer

*/
int check_for_patch(int iplume, struct plm *plume,
		    int ilayer,int apatch, struct lyr *layer, 
		    int nlon, int nlat,
		    GGRD_CPREC rmax,GGRD_CPREC *lon, 
		    GGRD_CPREC *lat,GGRD_CPREC z)
{
  int i,j,k,have_layer;
  if(ilayer <= 0){
    fprintf(stderr,"check_for_patch: should only use > 0\n");
    exit(-1);
  }
  /* check if this plume has any entries in the layer below ilayer */
  have_layer = -1;
  for(i=0;i<plume[iplume].nz;i++)
    if(plume[iplume].iz[i] == ilayer-1){
      have_layer = i;break;
    }
  if(have_layer != -1){
    /* 
       the layer below is in plume, check if the patch location is represented
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
    plume[iplume].iz = (int *)realloc(plume[iplume].iz,(plume[iplume].nz+1)*sizeof(int));
    plume[iplume].npixel = (int *)realloc(plume[iplume].npixel,(plume[iplume].nz+1)*sizeof(int));
    /*  */
    plume[iplume].pnr = (int *)realloc(plume[iplume].pnr,(plume[iplume].nz+1)*sizeof(int));
    plume[iplume].pname = (int **)realloc(plume[iplume].pname,(plume[iplume].nz+1)*sizeof(int *));
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
    /* didn't find this patch code, add to list */
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
	plume[iplume].pixels[luse][plume[iplume].npixel[luse]].c[0] = j;
	plume[iplume].pixels[luse][plume[iplume].npixel[luse]].c[1] = i;
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

/* 

compute the mean conduit properties for all plumes 

*/
void compute_mean_plume_properties(struct plm *plume, int nplume,GGRD_CPREC *w, GGRD_CPREC ws,
				   struct lyr *layer)
{
  int i,j,k;
  for(i=0;i < nplume;i++){	/* loop through plumes */
    /* make room for area computation */
    plume[i].area = (GGRD_CPREC *)realloc(plume[i].area,sizeof(GGRD_CPREC)*plume[i].nz);
    if(!plume[i].area){fprintf(stderr,"mem error\n");exit(-1);}
    plume[i].zmin = 1e20;plume[i].zmax = -1e20;
    for(j=0;j < plume[i].nz;j++){	/* loop through layers */
      plume[i].mean_npixel += plume[i].npixel[j];
      plume[i].area[j] = 0.0;
      for(k=0;k < plume[i].npixel[j];k++)
	plume[i].area[j] += w[plume[i].pixels[j][k].c[1]];
      plume[i].area[j] /= ws;	/* area of plume entries for this layer */
      plume[i].mean_area += plume[i].area[j];
      if(layer[plume[i].iz[j]].z > plume[i].zmax)
	plume[i].zmax = layer[plume[i].iz[j]].z;
      if(layer[plume[i].iz[j]].z < plume[i].zmin)
	plume[i].zmin = layer[plume[i].iz[j]].z;
    }
    /* mean number of pixel */
    plume[i].mean_npixel = (int)((float)plume[i].mean_npixel/(float)plume[i].nz+.5);
    plume[i].mean_area /= (float)plume[i].nz;
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
  (*plume+(*nplume))->nz = 0;
  (*plume+(*nplume))->iz = NULL;
  (*plume+(*nplume))->pixels = NULL;
  (*plume+(*nplume))->npixel = NULL;
  (*plume+(*nplume))->pname = NULL;
  (*plume+(*nplume))->pnr = NULL;
  (*plume+(*nplume))->area = NULL;
  (*plume+(*nplume))->mean_area = 0.0;
  (*plume+(*nplume))->mean_npixel = 0;
  if(!(*plume)){
    fprintf(stderr,"mem error plume %i\n",*nplume+1);
    exit(-1);
  }
  *nplume += 1;
}
/* in km */
GGRD_CPREC int_dist(int i1,int j1,int i2,int j2,
	       GGRD_CPREC *lon,GGRD_CPREC *lat,
	       GGRD_CPREC z)
{
  GGRD_CPREC lon1,lat1,lon2,lat2,r,tmp1,tmp2,tmp3;
  static GGRD_CPREC pif = GGRD_PI/180;
  r = 6371-fabs(z);		/* km */
  
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
