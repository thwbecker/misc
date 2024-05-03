/* 
*/
#include <math.h>
#include "gmt.h"

#define GMT_PRECISION float
float objective_func(float *, struct GRD_HEADER *, float *,int);
float fitf(float *, float, float, int);
void amoeba(float **, float *, int, float, float (*)(void), int *, struct GRD_HEADER *, float *,int);
float ran3(long *);
float amotry(float **, float *, float *, int, float (*)(void), int, float, struct GRD_HEADER *, float *,int);
float *vector(long, long);
void free_vector(float *, long, long);
float **matrix(long, long, long, long);
void free_matrix(float **, long, long, long, long);

  
#define NR_END 1
#define FREE_ARG char*

#define NPAR 4

int main(int argc, char **argv)
{
  struct GRD_HEADER header[1];
  GMT_PRECISION *val,eps,x,y,*gmt_val,max,min,xmax,ymax;
  GMT_LONG dummy[4]={0,0,0,0};
  /* 
   A = par[1]
   x_0 = par[2]
   y_0 = par[3]
   sigma = sqrt(par[4]/2)
 */
  GMT_PRECISION **p,y0[NPAR+1]={0,
				1, 0, 0, 200}; /* start values */
  GMT_PRECISION par[NPAR+1],objective0[NPAR+2],misfit;
  long seed = -1;
  int ndim,nrestart=2,i,j,k,nfunk;
  int fmode = 3;		/* 1: Gaussian (exp(-x**2) 
				   2: Cosine(x)^2
				   3: only Pos gaussian
				 */
  FILE *out;
  
  ran3(&seed);
  GMT_program=argv[0];

  GMT_begin (argc, argv);
  GMT_grd_init (header, argc, argv, FALSE);

  /* 
     default values 
  */
  if(argc < 3){
    fprintf(stderr,"usage:\n%s file.grd out.file\n\n",
	    argv[0]);
    exit(-1);
  }
  /* 
     
  read grid with F(lon,lat) 

  */

  if(GMT_read_grd_info (argv[1],header)== -1){
    fprintf(stderr,"%s: error opening %s for grdinfo\n",argv[0],argv[1]);
    exit(-1);
  }
  fprintf(stderr,"%s: read info OK\n",argv[0]);
  /* cartesian */

  //header->x_min;     
  //header->x_max;
  
  //header->y_min;
  //header->y_max;
  
  // header->x_inc;
  // header->y_inc;
  
  // header->nx;
  // header->ny;
  /* 
     
  read the grid

  */
  fprintf(stderr,"%s: reading from grd file %s...\n",argv[0],argv[1]);
  if((gmt_val=(GMT_PRECISION *)malloc(sizeof(GMT_PRECISION)*(header->nx+2)*(header->ny+2)))==NULL){
    fprintf(stderr,"%s: memerror, val too large: %i by %i\n",
	    argv[0],header->nx,header->ny);
    exit(-1);
  }
  GMT_read_grd (argv[1],header,gmt_val, 0,0,0,0,dummy,0);
  /* resort */
  if((val=(GMT_PRECISION *)malloc(sizeof(GMT_PRECISION)*(header->nx)*(header->ny)))==NULL){
    fprintf(stderr,"%s: memerror, val too large: %i by %i\n",
	    argv[0],header->nx,header->ny);
    exit(-1);
  }
  /* resort */
  max=-1e20;min=1e20;xmax=ymax=1e20;
  for(i=0;i < header->ny;i++){
    y = header->y_min + ((GMT_PRECISION)i)*header->y_inc;
    for(j=0;j < header->nx;j++){
      x = header->x_min + ((GMT_PRECISION)j)*header->x_inc;
      val[i*header->nx+j] = gmt_val[(header->ny-1-i)*header->nx+j];
      if(val[i*header->nx+j] < min)
	min = val[i*header->nx+j]; /*  */
      if(val[i*header->nx+j] > max){
	max = val[i*header->nx+j]; /*  */
	xmax = x; ymax = y;
      }
    }
  }
  free(gmt_val);
  fprintf(stderr,"%s: read grid: xrange:(%g %g %g) yrange (%g %g %g): min %g max %g at (%g,%g)\n",
	  argv[1],
	  header->x_min,header->x_inc,header->x_max,
	  header->y_min,header->y_inc,header->y_max,
	  min,max,xmax,ymax);
  
  y0[1] = max;y0[2] = xmax;y0[3] = ymax; /* starting guesses */

  /* invert */
  ndim = NPAR;
  eps = 0.1;			/* init offset */
  
  p = matrix(1,ndim+1,1,ndim);

  for(i=0;i<nrestart;i++){
    for(j=1;j<=ndim+1;j++){	/* init */
      for(k=1;k<=ndim;k++){
	if(i==0)		/* initial call */
	  par[k] = y0[k];
	else			/* restart */
	  par[k] = p[1][k] + (-.5+ran3(&seed));
	if((j != 1)&&(k==j-1))
	  par[k] += eps;
	p[j][k] = par[k];
      }
      objective0[j] = objective_func(par, header, val,fmode);
      if(i==0)
	fprintf(stderr,"%s: init: %12g %12g %12g %12g, misfit: %7.3e\n",
		argv[0],
		p[j][1], p[j][2], p[j][3], p[j][4],
		objective0[j]);
      
    }
    
    amoeba(p,objective0,ndim,1e-5,(float (*)(void))objective_func,&nfunk,header,val,fmode);
    misfit = objective_func(p[1], header, val,fmode);
    fprintf(stderr,"%s: solve %i: %12g %12g %12g %12g  after %5i evaluations, misfit: %7.3e\n",
	    argv[0],k+1,p[1][1],p[1][2],p[1][3],p[1][4],nfunk,misfit);
  }
  /* 

     output 
  
  */
  if(0){
  }
  fprintf(stdout,"%g %g %g %g %g\n",p[1][1],p[1][2],p[1][3],p[1][4],misfit);
  /* fit output */
  out = fopen(argv[2],"w");
  /* x y data fit data-fit */
  for(i=0,y=header->y_min;y <= header->y_max;i++,y+=header->y_inc)
    for(j=0,x=header->x_min;x <= header->x_max;x+=header->x_inc,j++)
      fprintf(out,"%g %g %g %g %g\n",x,y,val[i*header->nx+j],fitf(p[1],x,y,fmode),
	      val[i*header->nx+j]-fitf(p[1],x,y,fmode));
  fclose(out);
  
  /* free space */
  free(val);
  free_matrix(p,1,ndim+1,1,ndim);

  return 0;
}


/* objective function to compute misfit with 2D Gaussian at */
GMT_PRECISION objective_func(GMT_PRECISION *par, struct GRD_HEADER *header, GMT_PRECISION *val,int fmode)
{
  int i,j,n;
  GMT_PRECISION res,x,y,diff;
  res = 0;
  n=0;
  for(i=0;i<header->ny;i++){
    y = header->y_min + ((GMT_PRECISION)i)*header->y_inc; /* y loc */
    for(j=0;j<header->nx;j++){
      x = header->x_min + ((GMT_PRECISION)j)*header->x_inc; /* x loc */
      //fprintf(stderr,"%g %g %g %g\n",x,y,fitf(par,x,y,fmode),val[i*header->nx+j]);
      if(finite(val[i*header->nx+j])){
	diff = fitf(par,x,y,fmode) - val[i*header->nx+j];
	res += diff*diff;
	n++;
      }
    }
  }
  if(fmode == 3)
    if(par[1]<0)
      res += 1e20;
  return res/(GMT_PRECISION)n;
}


/* 
   fitting function

   mode = 1
   
   A = par[1]
   x_0 = par[2]
   y_0 = par[3]
   sigma = sqrt(par[4])


 */
GMT_PRECISION fitf(GMT_PRECISION *par, GMT_PRECISION x, GMT_PRECISION y,int fmode)
{
  double xd,yd,val,dist,tmp;
  xd = (double)x - (double)par[2];
  yd = (double)y - (double)par[3];
  
  dist = (xd*xd+yd*yd)/par[4];	/* distance**2/par[4] */
  if((fmode == 1)||(fmode == 3)){
    val = (double)par[1] * exp(-dist);
  }else if(fmode == 2){
    tmp = cos(dist);
    val = (double)par[1] *tmp*tmp;
  }else{
    fprintf(stderr,"fitf: mode error %i",fmode);
    exit(-1);
  }
  return (float)val;

}



#define NMAX 5000
#define GET_PSUM \
					for (j=1;j<=ndim;j++) {\
					for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
					psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void amoeba(GMT_PRECISION **p,GMT_PRECISION *y,int ndim,GMT_PRECISION ftol,GMT_PRECISION (*funk)(),
	    int *nfunk, struct GRD_HEADER *header, GMT_PRECISION *val,int fmode)
{
  int i,ihi,ilo,inhi,j,mpts=ndim+1;
  GMT_PRECISION rtol,sum,swap,ysave,ytry,*psum;

	psum=vector(1,ndim);
	*nfunk=0;
	GET_PSUM
	for (;;) {
		ilo=1;
		ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
		for (i=1;i<=mpts;i++) {
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi) inhi=i;
		}
		rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
		if (rtol < ftol) {
			SWAP(y[1],y[ilo])
			for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
			break;
		}
		if (*nfunk >= NMAX) {
		  fprintf(stderr,"NMAX exceeded");
		  exit(-1);
		}
		*nfunk += 2;
		ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0,header,val,fmode);
		if (ytry <= y[ilo])
		  ytry=amotry(p,y,psum,ndim,funk,ihi,2.0,header,val,fmode);
		else if (ytry >= y[inhi]) {
			ysave=y[ihi];
			ytry=amotry(p,y,psum,ndim,funk,ihi,0.5,header,val,fmode);
			if (ytry >= ysave) {
				for (i=1;i<=mpts;i++) {
					if (i != ilo) {
						for (j=1;j<=ndim;j++)
							p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
						y[i]=(*funk)(psum,header,val,fmode);
					}
				}
				*nfunk += ndim;
				GET_PSUM
			}
		} else --(*nfunk);
	}
	free_vector(psum,1,ndim);
}
#undef SWAP
#undef GET_PSUM
#undef NMAX

GMT_PRECISION amotry(GMT_PRECISION **p,GMT_PRECISION *y,GMT_PRECISION *psum,
		     int ndim,GMT_PRECISION (*funk)(),int ihi,GMT_PRECISION fac,
		     struct GRD_HEADER *header, GMT_PRECISION *val,int fmode)
{
	int j;
	GMT_PRECISION fac1,fac2,ytry,*ptry;

	ptry=vector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(ptry,header,val,fmode);
	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=1;j<=ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free_vector(ptry,1,ndim);
	return ytry;
}


GMT_PRECISION *vector(long nl,long nh)
/* allocate a GMT_PRECISION vector with subscript range v[nl..nh] */
{
	GMT_PRECISION *v;

	v=(GMT_PRECISION *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(GMT_PRECISION)));
	if (!v) {
	  fprintf(stderr,"allocation failure in vector()");
	  exit(-1);
	}
	return v-nl+NR_END;
}

void free_vector(GMT_PRECISION *v,long nl,long nh)

/* free a GMT_PRECISION vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
GMT_PRECISION **matrix(long nrl,long nrh,long ncl,long nch)
/* allocate a GMT_PRECISION matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	GMT_PRECISION **m;

	/* allocate pointers to rows */
	m=(GMT_PRECISION **) malloc((unsigned int)((nrow+NR_END)*sizeof(GMT_PRECISION*)));
	if (!m) {
	  fprintf(stderr,"allocation failure 1 in matrix()");
	  exit(-1);
	}
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(GMT_PRECISION *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(GMT_PRECISION)));
	if (!m[nrl]) {
	  fprintf(stderr,"allocation failure 2 in matrix()");
	  exit(-1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}
void free_matrix(GMT_PRECISION **m,long nrl,long nrh,long ncl,long nch)
/* free a GMT_PRECISION matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(long *idum)
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

