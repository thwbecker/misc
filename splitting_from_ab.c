#include "splitting_from_ab.h"
/* 

read a set of a, b coefficient grd files from a 3D azimuthal
anisotropy model and predict splitting

*/

int main(int argc, char **argv)
{
  char dir_name[HC_CHAR_LENGTH];
  double lon,lat,*dweight,*weight,phi,dt,zrange,zmin,gsum[2];
  int i,nweight,wmode,n,nweight_0;
  struct mod model[1];
  
  nweight_0 = 150;		/* subdivisions for integrals */

  wmode = 1;			/* weight mode */
  switch(argc){
  case 2:
    sprintf(dir_name,"%s",argv[1]);
    break;
  case 3:
    sprintf(dir_name,"%s",argv[1]);
    sscanf(argv[2],"%i",&wmode);
    break;
     
  default:
    fprintf(stderr,"%s: usage\n%s dir_name [wmode, %i]\n\nreads a, b grd files from directory dir_name\n",
	    argv[0],argv[0],wmode);
    exit(-1);
    break;
  }
  /* init PREM */
  prem_read_model(PREM_MODEL_FILE,model->prem,TRUE);

  /* init grids */
  read_ab(model,dir_name);
  
  /*  */
  switch(wmode){
  case 0:			/* original */
    /* assign weights */
    nweight = model->nz;
    hc_dvecalloc(&weight,nweight,argv[0]);
    hc_dvecalloc(&dweight,nweight,argv[0]);
    for(i=0;i < nweight;i++){
      weight[i] = 1.0;
      dweight[i] = -model->a[0].z[i];
    }
    break;
  case 1:			/* single layer at 150 km*/
    nweight = 2;zrange = 15;zmin = 150+zrange/2;
    hc_dvecalloc(&weight,nweight,argv[0]);
    hc_dvecalloc(&dweight,nweight,argv[0]);
    for(i=0;i < nweight;i++){
      weight[i] = 1.0;
      dweight[i] = zmin - ((double)i/(double)(nweight-1))*zrange;
    }
    break;
  case 2:			/* between 410 and 75 */
    nweight = nweight_0;zrange = 335;zmin = 410;
    hc_dvecalloc(&weight,nweight,argv[0]);
    hc_dvecalloc(&dweight,nweight,argv[0]);
    for(i=0;i < nweight;i++){
      weight[i] = 1.0;
      dweight[i] = zmin - ((double)i/(double)(nweight-1))*zrange;
    }
    break;
  case 3:			/* between 660 and 75 */
    nweight = nweight_0;zrange = 585;zmin = 660;
    hc_dvecalloc(&weight,nweight,argv[0]);
    hc_dvecalloc(&dweight,nweight,argv[0]);
    for(i=0;i < nweight;i++){
      weight[i] = 1.0;
      dweight[i] = zmin - ((double)i/(double)(nweight-1))*zrange;
    }
    break;
  case 4:			/* between 410 and 10 */
    nweight = nweight_0;zrange = 400;zmin = 410;
    hc_dvecalloc(&weight,nweight,argv[0]);
    hc_dvecalloc(&dweight,nweight,argv[0]);
    for(i=0;i < nweight;i++){
      weight[i] = 1.0;
      dweight[i] = zmin - ((double)i/(double)(nweight-1))*zrange;
    }
    break;
  case 5:			/* between 660 and 25 */
    nweight = nweight_0;zrange = 635;zmin = 660;
    hc_dvecalloc(&weight,nweight,argv[0]);
    hc_dvecalloc(&dweight,nweight,argv[0]);
    for(i=0;i < nweight;i++){
      weight[i] = 1.0;
      dweight[i] = zmin - ((double)i/(double)(nweight-1))*zrange;
    }
    break;
  case 6:			/* between 250 and 25 */
    nweight = nweight_0;zrange = 225;zmin = 250;
    hc_dvecalloc(&weight,nweight,argv[0]);
    hc_dvecalloc(&dweight,nweight,argv[0]);
    for(i=0;i < nweight;i++){
      weight[i] = 1.0;
      dweight[i] = zmin - ((double)i/(double)(nweight-1))*zrange;
    }
    break;
  case 7:			/* between 250 and 75 */
    nweight = nweight_0;zrange = 175;zmin = 250;
    hc_dvecalloc(&weight,nweight,argv[0]);
    hc_dvecalloc(&dweight,nweight,argv[0]);
    for(i=0;i < nweight;i++){
      weight[i] = 1.0;
      dweight[i] = zmin - ((double)i/(double)(nweight-1))*zrange;
    }
    break;
  case 8:			/* between 110 and 25 */
    nweight = nweight_0;zrange = 85;zmin = 110;
    hc_dvecalloc(&weight,nweight,argv[0]);
    hc_dvecalloc(&dweight,nweight,argv[0]);
    for(i=0;i < nweight;i++){
      weight[i] = 1.0;
      dweight[i] = zmin - ((double)i/(double)(nweight-1))*zrange;
    }
    break;
  default:
    fprintf(stderr,"%s: wmode %i undefined\n",argv[0],wmode);
    exit(-1);
    break;
  }
  n=0;
  fprintf(stderr,"%s: reading lon lat from stdin, writing lon lat phi[deg] dt[s] dsum[0] dsum[1]\n",
	  argv[0]);
  while(fscanf(stdin,"%lf %lf",&lon,&lat)==2){
    compute_sks(model,lon,lat,nweight,dweight,weight,&phi,&dt,gsum,
		(n==0));
    fprintf(stdout,"%11g %11g %11g %11g %11g %11g\n",
	    lon,lat,RAD2DEG(phi),dt,gsum[0],gsum[1]);
    n++;
  }
  fprintf(stderr,"%s: read %i locations\n",argv[0],n);
  return 0;
}

/* 

   compute SKS fast angle phi [rad] and delay time dt [sec] given
   geographic location lon/lat [deg] and weights specified on nweight
   layers, with depth dweight [km, positive]


 */
void compute_sks(struct mod *model, double lon, double lat,
		 int nweight, double *dweight, double *weight,
		 double *phi, double *dt,double *gsum,
		 hc_boolean show_sum)
{
  int i,j,k;
  double z,dummy[2],vsv,x_r,x_theta,x_phi,a[2],dz;
  *dt = 0;
  *phi = 0;

  x_theta = LAT2THETA(lat);
  x_phi =   LON2PHI(lon);
  gsum[0] = gsum[1] = 0.0;
  if(nweight < 2){
    fprintf(stderr,"compute_sks: need at least two depth levels (defining one layer)\n");
    exit(-1);
  }
  for(i=0;i < nweight-1;i++){
    dz = -(dweight[i+1]-dweight[i]);     /* layer thickness */
    /* mid layer depth */
    z = -(dweight[i]+dweight[i+1])/2;		/* z counts negative downward */

    if(z >= 0){
      fprintf(stderr,"compute_sks: error: layer %i, z: %g depth: %g\n",
	      i+1,z,dweight[i]);
      exit(-1);
    }
    /* obtain vsv in km/s at mid layer depths */
    prem_get_values(dummy,dummy,&vsv,dummy,dummy,dummy,dummy,
		    dummy,dummy,dummy,(HC_RE_KM + z)*1e3, 
		    model->prem);vsv /= 1000.;
    /* get A/B at this layer and above */
    a[0] = a[1] = 0.0;
    for(j=0;j < 2;j++){
      x_r = HC_ND_RADIUS(dweight[i+j]);
      for(k=0;k < 2;k++){
	ggrd_grdtrack_interpolate_rtp(x_r,x_theta,x_phi,&(model->a[k]),
				      (dummy+k),TRUE,FALSE,6371.);
	dummy[k] /= 50;		/* assume anomalies are in % */
	a[k] += dummy[k];		/* G/L = 2A */
      }
      //fprintf(stderr,"d: %g a: %g b: %g\n",dweight[i+j],dummy[0],dummy[1]);
    }
    for(j=0;j<2;j++){		
      a[j] /= 2;		/* get mean, assuming change is linear */
      gsum[j] += dz * a[j]/vsv;	/* add up */
    }
    if(show_sum)
      fprintf(stderr,"d1: %6.1f d2: %6.1f mid: %6.2f km dz: %7.3f weight: %6.4f vsv: %7.3f [km/s] a: %8.5f b: %8.5f\n",
	      dweight[i],dweight[i+1],-z,dz,weight[i],vsv,a[0],a[1]);
  }
  *dt  = hypot(gsum[0],gsum[1]);
  *phi = 0.5*atan2(gsum[1],gsum[0]);
  if(*phi < 0)
    *phi += HC_TWOPI;
  if(*phi > HC_PI)
    *phi -= HC_PI;
}
/* 

read global a,b gxrids in % wrt to PREM

*/
void read_ab(struct mod *model, char *dir_name)
{
  char grd_name[HC_CHAR_LENGTH],dfile_name[HC_CHAR_LENGTH];
  ggrd_boolean verbose = TRUE;
  sprintf(dfile_name,"%s/depths.dat",dir_name);

  sprintf(grd_name,"%s/a",dir_name);
  if(ggrd_grdtrack_init_general(TRUE,grd_name,dfile_name,"-fg",
				(model->a+0),verbose,TRUE,FALSE)){
    fprintf(stderr,"read error\n");
    exit(-1);
  }

  sprintf(grd_name,"%s/b",dir_name);
  if(ggrd_grdtrack_init_general(TRUE,grd_name,dfile_name,"-fg",
				(model->a+1),verbose,TRUE,FALSE)){
    fprintf(stderr,"read error\n");
    exit(-1);
  }

  /* depth counted negative downward */
  model->nz = model->a[0].nz;
  model->top    = model->a[0].z[model->a[0].nz-1];
  model->bottom = model->a[0].z[0];
  fprintf(stderr,"read_ab: read model, %i layers, bottom: %g top: %g\n",
	  model->nz,model->bottom, model->top);
}
