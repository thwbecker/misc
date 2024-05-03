/* 

solve poisson's equation for a grd file which holds F(lon,lat)

either on spherical, or Cartesian geometry

if gbounded, can set either fixed, zero value or fixed, zero gradient
as boundary condition, or free-slip like (assuming field is for
velocities)

uses FISHPAK


$Id: grd2poisson.c,v 1.3 2008/05/16 14:33:47 twb Exp twb $

*/
#ifdef USE_GMT4
#include "gmt.h"
#else
#include "gmt_dev.h"
//#include "gmt_private.h"
#endif

#define GMT_PRECISION float
#include "fishpak.h"
/* trig stuff */
#define PI 3.14159265358979324
#define TWOPI 6.28318530717958647
#define PIOVERONEEIGHTY 0.0174532925199433
#define ONEEIGHTYOVERPI  57.2957795130823
#define THETA2LAT(x) ( (90.0 - (x)*ONEEIGHTYOVERPI) )
#define LAT2THETA(x) ( (90.0 - (x))*PIOVERONEEIGHTY )
#define PHI2LON(x) ( (x) * ONEEIGHTYOVERPI )
#define LON2PHI(x) ( (x) * PIOVERONEEIGHTY )


int main(int argc, char **argv)
{

  int sp=0,np=0,phiper=0,thetaper=0,n,m,
    phibc,thetabc,spherical,nphi,ntheta,i,j,mbcond,nbcond,ierror,idimf;
  GMT_PRECISION *func,theta_min,theta_max,phi_min,phi_max,dtheta,dphi;
  char ofile[2000];
  FISH_PREC a,b,*bda,*bdb,c,d,*bdc,*bdd,lambda,
    *w,pertrb,*f;

  char *cdummy[1];
  FILE *in;
#ifdef USE_GMT4
  struct GRD_HEADER header[1];
  GMT_LONG dummy[4]={0,0,0,0};
  GMT_program=argv[0];
  GMT_begin (argc, argv);
  GMT_grd_init (header, argc, argv, FALSE);
#else
  void *API;                        /* The API control structure */
  struct GMT_GRID *G = NULL;        /* Structure to hold output grid */
  struct GMT_GRID_HEADER *header = NULL;
  int ij,pixelreg;
  double wesn[6];
  double inc[2];
  int pad;
  uint64_t par[2];
  API = GMT_Create_Session (argv[0], 2U, 0, NULL);
#endif
  /* 
     default values 
  */
  spherical = 1;
  phibc = -1;
  phiper = thetaper = 0;


  if(argc < 2){
    fprintf(stderr,"usage:\n%s file.grd [bcx, %i] [bcy, bcx] [s, %i]\ncompute solution to Poisson's eq\n\n",
	    argv[0],phibc,spherical);
    fprintf(stderr,"bc -1: gradients specified in files a.dat b.dat c.dat d.dat (bottom up, left right)\n");
    fprintf(stderr,"    0: zero gradients on boundary\n");
    fprintf(stderr,"    1: zero value     on boundary\n");
    fprintf(stderr,"    2: periodic\n");
    
    fprintf(stderr,"s   1: spherical coordinates\n");
    fprintf(stderr,"    0: Cartesian coordinates\n");
    exit(-1);
  }
  /* which boundary condition */
  if(argc > 2)
    sscanf(argv[2],"%i",&phibc);
  /* use same boundary conditions in both  */
  thetabc = phibc;		
  if(argc > 3)
    sscanf(argv[3],"%i",&thetabc);
  /* 
     spherical system or cartesian
  */
  if(argc > 4)
    sscanf(argv[4],"%i",&spherical);
  if(spherical)
    fprintf(stderr,"%s: Spherical geometry\n",argv[0]);
  else 
    fprintf(stderr,"%s: Cartesian geometry\n",argv[0]);
  /* 
     
  read grid with F(lon,lat) 

  */
#ifdef USE_GMT4
  if(GMT_read_grd_info (argv[1],header)== -1){
    fprintf(stderr,"%s: error opening %s for grdinfo\n",argv[0],argv[1]);
    exit(-1);
  }
  fprintf(stderr,"%s: read info OK\n",argv[0]);
  if(!spherical){
    /* cartesian */
    sp = np = 0; /* no periodicity by default */
    if(phibc == 2)
      phiper = 1;
    if(thetabc == 2)
      thetaper = 1;
    
    phi_min = header->x_min;     
    phi_max = header->x_max;

    theta_min = header->y_min;
    theta_max = header->y_max;

    dphi   = header->x_inc;
    dtheta = header->y_inc;

    nphi   = header->nx;
    ntheta = header->ny;

  }else{
    thetaper = 0;		/* cannot be periodic in theta */
    if(fabs(header->y_min+90)<1e-5){
      sp = 1;			/* south pole included */
      theta_max = pimach();
    }else{
      theta_max = LAT2THETA(header->y_min);
    }
    if(fabs(header->y_max-90)<1e-5){
      np = 1;			/* north pole included */
      theta_min = 0.0;
    }else{
      theta_min=LAT2THETA(header->y_max);
    }
    phi_min = LON2PHI(header->x_min);
    phi_max = LON2PHI(header->x_max);

    if(phibc == 2)
      phiper = 1;
    if(thetabc == 2){
      fprintf(stderr,"%s: theta periodic doesn't make sense for spherical\n",argv[0]);
      exit(-1);
    }
    
    dphi   = LON2PHI(header->x_inc);
    dtheta = LON2PHI(header->y_inc);

    nphi =   header->nx;
    ntheta = header->ny;
  }
  /* 
     
  read the grid

  */
  fprintf(stderr,"%s: reading from grd file %s...\n",argv[0],argv[1]);
  if((func=(GMT_PRECISION *)malloc(sizeof(GMT_PRECISION)*(ntheta+2)*(nphi+2)))==NULL){
    fprintf(stderr,"%s: memerror, func too large: %i by %i\n",
	    argv[0],nphi,ntheta);
    exit(-1);
  }
  GMT_read_grd (argv[1],header,func, 0,0,0,0,dummy,0);
  fprintf(stderr,"%s: read grid OK\n",argv[1]);
#else
  /* GMT 6 version  */
  if((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE,
			 GMT_IS_SURFACE, GMT_CONTAINER_ONLY,
			 NULL, argv[1], NULL))==NULL)
    return (-1);
  header = G->header;
  pixelreg=(header->registration == GMT_GRID_PIXEL_REG)?1:0;
    
  phi_min=LON2PHI(header->wesn[XLO]+(pixelreg?header->inc[GMT_X]/2.0:0.0));	
  phi_max=LON2PHI(header->wesn[XHI]-(pixelreg?header->inc[GMT_X]/2.0:0.0));
  theta_min=LAT2THETA(header->wesn[YLO]+(pixelreg?header->inc[GMT_Y]/2.0:0.0));	
  theta_max=LAT2THETA(header->wesn[YHI]-(pixelreg?header->inc[GMT_Y]/2.0:0.0));
  dphi=  LON2PHI(header->inc[GMT_X]);
  dtheta=LON2PHI(header->inc[GMT_Y]);
  // for grid line registered nx_i = (x_i^max-x_i^min)/dx_i + 1
  // for pixel reg            nx_i = (x_i^max-x_i^min)/dx_i
  nphi=header->n_columns;
  ntheta=header->n_rows;

  if(!spherical){
    /* cartesian */
    sp = np = 0; /* no periodicity by default */
    if(phibc == 2)
      phiper = 1;
    if(thetabc == 2)
      thetaper = 1;
  }else{
    thetaper = 0;		/* cannot be periodic in theta */
    if(fabs(theta_max+90)<1e-5){
      sp = 1;			/* south pole included */
      theta_max = pimach();
    }
    if(fabs(theta_min-90)<1e-5){
      np = 1;			/* north pole included */
      theta_min = 0.0;
    }
    if(phibc == 2)
      phiper = 1;
    if(thetabc == 2){
      fprintf(stderr,"%s: theta periodic doesn't make sense for spherical\n",argv[0]);
      exit(-1);
    }
    
  }
  /* 
     
  read the grid

  */
  fprintf(stderr,"%s: reading from grd file %s...\n",argv[0],argv[1]);
  if((func=(GMT_PRECISION *)malloc(sizeof(GMT_PRECISION)*(ntheta+2)*(nphi+2)))==NULL){
    fprintf(stderr,"%s: memerror, func too large: %i by %i\n",
	    argv[0],nphi,ntheta);
    exit(-1);
  }
   if((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE,
			  GMT_IS_SURFACE, GMT_DATA_ONLY, NULL, argv[1], G))==NULL)
	return (-1);
  fprintf(stderr,"%s: read grid OK\n",argv[1]);
#endif

  /*
    
  prepare call to fishpak

  */
  if(!spherical){
    /* 
       cartesian 

    */

    /* phi direction */
    a = phi_min;
    b = phi_max;

    m = nphi-1;
    n = ntheta-1;
    /* 
       
    phi BC
    
    */
    if(phiper)
      mbcond = 0;		/* x periodic */
    else{
      switch(phibc){
      case 1:
	mbcond = 1;			/* all BC values fixed at zero */
	break;
      case  0:
      case -1:
	mbcond = 3;			/* derivatives set to zero or value */
	break;
      default:
	fprintf(stderr,"%s: bc error, %i undefined\n",argv[0],phibc);
	exit(-1);
	break;
      }
    }
    /* init */
    bda = (FISH_PREC *)calloc(ntheta,sizeof(FISH_PREC));
    bdb = (FISH_PREC *)calloc(ntheta,sizeof(FISH_PREC));
    /* 
       phi direction 
    */
    c = theta_min;
    d = theta_max;
    if(thetaper)
      nbcond = 0;		/* y periodic */
    else{
      switch(thetabc){
      case 1:
	nbcond = 1;			/* all BC values fixed at zero */
	break;
      case  0:
      case -1:
	nbcond = 3;			/* derivatives set to zero or value */
	break;
      default:
	fprintf(stderr,"bc error\n");
	exit(-1);
	break;
      }
    }
    /* init as zeros */
    bdc = (FISH_PREC *)calloc(nphi,sizeof(FISH_PREC));
    bdd = (FISH_PREC *)calloc(nphi,sizeof(FISH_PREC));
    lambda = 0.0;
    /*  */
    f = (FISH_PREC *)calloc(nphi*ntheta,sizeof(FISH_PREC));
    w = (FISH_PREC *)calloc(20*ntheta+4*nphi+ntheta*nphi,
			  sizeof(FISH_PREC));
    if(!bda || !bdb || !bdc || !bdd || !f || !w){
      fprintf(stderr,"%s: mem error\n",argv[0]);
      exit(-1);
    }
    /* 
       gmt sorting is
       x_i, y_i --> x((ny-j-1)*nx+i) 
    */
   
    /* 
       resort for fortran?
       
       F(I,J) = F(X(I),Y(J)) = j*nx + i
       
    */
    for(i=0;i<ntheta;i++)
      for(j=0;j<nphi;j++)
	f[i*nphi+j] = func[(ntheta-1-i)*nphi+j];
    idimf = nphi;
    if(phibc == -1){
      /* read gradient values from file */
      in = fopen("a.dat","r");if(!in){fprintf(stderr,"%s: could not open a.dat\n",argv[0]);exit(-1);}
      for(i=0;i<ntheta;i++)
	if(fscanf(in,"%f",(bda+i))!=1){fprintf(stderr,"%s: read error a.dat\n",argv[0]);exit(-1);}
      fclose(in);
      in = fopen("b.dat","r");if(!in){fprintf(stderr,"%s: could not open b.dat\n",argv[0]);exit(-1);}
      for(i=0;i<ntheta;i++)
	if(fscanf(in,"%f",(bdb+i))!=1){fprintf(stderr,"%s: read error b.dat\n",argv[0]);exit(-1);}
      fclose(in);
      /*  */
      in = fopen("c.dat","r");if(!in){fprintf(stderr,"%s: could not open c.dat\n",argv[0]);exit(-1);}
      for(i=0;i<nphi;i++)
	if(fscanf(in,"%f",(bdc+i))!=1){fprintf(stderr,"%s: read error c.dat\n",argv[0]);exit(-1);}
      fclose(in);
      in = fopen("d.dat","r");if(!in){fprintf(stderr,"%s: could not open d.dat\n",argv[0]);exit(-1);}
      for(i=0;i<nphi;i++)
	if(fscanf(in,"%f",(bdd+i))!=1){fprintf(stderr,"%s: read error d.dat\n",argv[0]);exit(-1);}
      fclose(in);
      fprintf(stderr,"%s: read all gradients in OK\n",argv[0]);


    }
    /*  */
    fprintf(stderr,"%s: a: %g b: %g c: %g d: %g\n",argv[0],a,b,c,d);
    fprintf(stderr,"%s: xbcond: %i ybcond: %i nx: %i ny: %i\n",argv[0],
	    mbcond,nbcond,nphi,ntheta);
    /* 
       
    solve cartesian
    
    */
    hwscrt(&a,&b,&m,&mbcond,bda,bdb,&c,&d,&n,&nbcond,bdc,bdd,
	   &lambda,f,&idimf,&pertrb,&ierror,w);
    if(ierror != 0){
      fprintf(stderr,"%s: hwscrt returned error code %i\n",
	      argv[0],ierror);
      exit(-1);
    }
    /* resort y axes */
    for(i=0;i<ntheta;i++)
      for(j=0;j<nphi;j++)
	func[(ntheta-1-i)*nphi+j] = f[i*nphi+j];
    /* done cartesian */

    if(0){
      for(i=0;i<ntheta;i++)
	if(mbcond == 3)
	  fprintf(stderr,"du/dx: y: %11g x = %11g: %11g <-- %11g x = %11g: %11g <-- %11g\n",
		  theta_max-dtheta*i,
		  phi_min,(func[i*nphi+1]-func[i*nphi])/dphi,bda[i],
		  phi_max,(func[i*nphi+nphi-1]-func[i*nphi+nphi-2])/dphi,bdb[i]);
	else
	  fprintf(stderr,"u: y: %11g x = %11g: %11g  x = %11g: %11g\n",
		  theta_max-dtheta*i,
		  phi_min,func[i*nphi],phi_max,func[i*nphi+nphi-1]);


      fprintf(stderr,"\n");
      for(i=0;i<nphi;i++)
	if(nbcond == 3)
	  fprintf(stderr,"du/dy: x: %11g y = %11g: %11g <-- %11g y = %11g: %11g <-- %11g\n",
		  phi_min+dphi*i,
		  theta_min,(func[(ntheta-1)*nphi+i]-func[(ntheta-2)*nphi+i])/dtheta,bdc[i],
		  theta_max,(func[                i]-func[nphi+i])/dtheta,bdd[i]);
	else
	  fprintf(stderr,"u: x: %11g y = %11g: %11g y = %11g: %11g\n",
		  phi_min+dphi*i,
		  theta_min,func[(ntheta-1)*nphi+i],theta_max,func[i]);
      
    }
    
    
  }else{
    /* 

    SPHERICAL SOLUTION

    */


    /* theta direction */
    a = theta_min;
    b = theta_max;
    /* 
       
    theta BC
    
    */
    if(np && sp)			/* both north and south pole */
      mbcond = 9;
    else if(np){			/* north pole */
      switch(thetabc){
      case 1:
	mbcond = 5;			/* solution */
	break;
      case 0:
	mbcond = 6;			/* derivative */
	break;
      default:
	fprintf(stderr,"%s: spherical thetabc %i not implemented\n",argv[0],thetabc);
	break;
      }
    }else if(sp){			/* south pole */
      switch(thetabc){
      case 1:
	mbcond = 7;			/* solution */
	break;
      case 0:
	mbcond = 8;			/* derivative */
	break;
      default:
	fprintf(stderr,"%s: spherical thetabc %i not implemented\n",argv[0],thetabc);
	break;
      }
    }else{			/* regular theta boundary 
				   condition */
     switch(thetabc){
      case 1:
	mbcond = 1;			/* solution */
	break;
      case 0:
	mbcond = 3;			/* derivative */
	break;
      default:
	fprintf(stderr,"%s: spherical thetabc %i not implemented\n",argv[0],thetabc);
	break;
      }
    }
    bda = (FISH_PREC *)calloc(ntheta,sizeof(FISH_PREC));
    bdb = (FISH_PREC *)calloc(ntheta,sizeof(FISH_PREC));
    /* 
       phi direction 
    */
    c = phi_min;
    d = phi_max;
    if(phiper){/* periodic */
      nbcond = 0;			/* phi BC */
    }else{
      /* gap in global phi coverage  */
      switch(phibc){
      case 1:
	nbcond = 1;			/* solution */
	break;
      case 0:
	nbcond = 3;			/* derivative */
	break;
      default:
	fprintf(stderr,"%s: spherical thetabc %i not implemented\n",argv[0],thetabc);
	break;
      }
    }
    bdc = (FISH_PREC *)calloc(nphi,sizeof(FISH_PREC));
    bdd = (FISH_PREC *)calloc(nphi,sizeof(FISH_PREC));
    lambda = 0.0;
    /*  */
    f = (FISH_PREC *)calloc(nphi*ntheta,sizeof(FISH_PREC));
    w = (FISH_PREC *)calloc(13*ntheta+4*nphi+ntheta*nphi,
			  sizeof(FISH_PREC));
    if(!bda || !bdb || !bdc || !bdd || !f || !w){
      fprintf(stderr,"%s: mem error\n",argv[0]);
      exit(-1);
    }
    /* 
       resort for fortran:
       
       F(I,J) = F(THETA(I),PHI(J))
       
    */
    for(i=0;i<ntheta;i++)
      for(j=0;j<nphi;j++)
	f[j*ntheta+i] = func[i*nphi+j];
    idimf = ntheta;
    /*  */
    fprintf(stderr,"%s: NP: %i SP: %i glon: %i a: %g b: %g c: %g d: %g\n",
	    argv[0],np,sp,phiper,a,b,c,d);
    fprintf(stderr,"%s: tbcond: %i pbcond: %i\n",argv[0],
	    mbcond,nbcond);
    /* 
       
    solve
    
    */
    
    /* spherical */
    hstssp(&a,&b,&ntheta,&mbcond,bda,bdb,
	   &c,&d,&nphi,  &nbcond,bdc,bdd,
	   &lambda,f,&idimf,&pertrb,&ierror,w);
    if(ierror != 0){
      fprintf(stderr,"%s: hstcrt returned error code %i\n",
	      argv[0],ierror);
      exit(-1);
    }
    /* resort */
    for(i=0;i<ntheta;i++)
      for(j=0;j<nphi;j++)
	func[i*nphi+j] = f[j*ntheta+i];
  
  }
  /* 

     output 
  
  */
  sprintf(ofile,"%s.p",argv[1]);
  fprintf(stderr,"%s: %g - %g - %g, %g - %g - %g, %i %i\n",
	  argv[0],phi_min,dphi,phi_max,theta_min,dtheta,theta_max,
	  nphi,ntheta);
#ifdef USE_GMT4
  GMT_write_grd (ofile,header,func,0,0,0,0, dummy,0);
#else
  pad = -1;
  inc[GMT_Y] = dphi;
  inc[GMT_X] = dtheta;
  par[0] = nphi;par[1] = ntheta;
  wesn[0] = phi_min;wesn[1] = phi_max;
  wesn[2] = theta_min; wesn[3] = theta_max;
  wesn[4]=wesn[5]=0;
  
  G = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE,GMT_CONTAINER_AND_DATA,
		       par,wesn,inc,GMT_GRID_NODE_REG,pad,NULL);
  if(!G){
    fprintf(stderr,"my_gmt_write_grd: error creating grid for GMT > 4 mode\n");
    exit(-1);
  }
  for(j=0;j < nphi;j++){
    for(i=0;i < ntheta;i++){
      ij = gmt_M_ijp (G->header, i, j);
      G->data[ij] = func[i*nphi+j];
    }
  }
  GMT_Write_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
		  GMT_CONTAINER_AND_DATA, NULL, ofile, G);
  GMT_Destroy_Data(API,G->data);
#endif
  fprintf(stderr,"%s: written to %s\n",argv[0],ofile);
  /* free space */
  free(bda);free(bdb);free(func);
  free(f);free(bdc);free(bdd),free(w);

#ifdef USE_GMT4
    GMT_end (argc, argv);
#else
    GMT_Destroy_Session (API);
#endif

    return 0;
}
