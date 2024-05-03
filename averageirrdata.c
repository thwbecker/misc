#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*


  average x y data in several files given as arguments 
  to the program like

  averageirrdata file1 file2 file3


  will determine some dx spacing between the minimum and 
  maximum x values of all files and then interpolate 
  for each file to obtain the mean and standard deviation 
  
  output is

  xi yi yis

  the x value, mean y value, and standard deviation around yi

  for all nx

  for compilation, set INTERPOLATION to 

  0: uses spline interpolation
  1: uses linear interpolation
  2: uses quadratic interpolation
  3: uses cubic interpolation 
  and so on


  $Id: averageirrdata.c,v 1.3 2005/02/23 03:04:21 becker Exp becker $

*/

#define AVG_ID_INTERPOLATION 0

void spline(double *,double *,int ,double ,double ,double *);
void splint(double *,double *,double *,int ,double ,double *);
void polint(double *,double *,int ,double ,double *,
	    double *);
void locate(double *,int ,double ,int *);
static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))


int main(int argc, char **argv)
{
  double *x,*y,xmin,xmax,dx,*xi,*yi,*yis,ytmp,*y2,tx,dy;
  int i,j,nfile,*n,order,nrp,nrpf,nx,nmax,index;
  FILE *in;
  order = AVG_ID_INTERPOLATION + 1;
  /*

    read in all x y pairs

  */
  n=(int *)malloc(sizeof(int));/* number of points per file n[0] is zero 
				    n[1] the number in the first file, and so on */
  x=(double *)malloc(sizeof(double));
  y=(double *)malloc(sizeof(double));
  nfile=0;
  nrp=0;
  nmax=0;
  if(argc < 2){
    fprintf(stderr,"%s file1 file2 file3 ...\n",argv[0]);
    fprintf(stderr,"outputs a mean of interpolations of the x-y pairs in file_n\n");
    fprintf(stderr,"output format is\n\txi yi yis\nwhere yi is the mean and yis the std\n");
    fprintf(stderr,"at each xi interpolated on nx values. nx is max(n(file))*2\n");
    exit(-1);
  }
  for(i=1;i<argc;i++){
    in=fopen(argv[i],"r");
    if(!in){
      fprintf(stderr,"%s: error, file %i: \"%s\" not found\n",
	      argv[0],i,argv[i]);
      exit(-1);
    }
    // read in x y values for this file
    nrpf=0;
    while(fscanf(in,"%lf %lf",(x+nrp),(y+nrp))==2){
      nrpf++;
      nrp++;
      x=(double *)realloc(x,(nrp+1)*sizeof(double));
      y=(double *)realloc(y,(nrp+1)*sizeof(double));
      if(!x || !y){
	fprintf(stderr,"%s: memory error: nrp: %i\n",argv[0],nrp);
	exit(-1);
      }
    }
    fclose(in);
    n=(int *)realloc(n,(nfile+1)*sizeof(int));
    n[nfile] = nrpf;
    nfile++;
    if(nrpf > nmax)
      nmax=nrpf;
    fprintf(stderr,"%s: read in %5i x y pairs from file %3i (name: %20s), total points now: %6i\n",
	    argv[0],n[nfile-1],nfile,argv[i],nrp);
  }
  /*


   determine number of interpolation points from nmax


  */
#if (AVG_ID_INTERPOLATION == 0)// spline
  nx = nmax*4;
#elif (AVG_ID_INTERPOLATION > 0)// polynomial
  nx = nmax*4;
#endif
  if(nx<4)
    nx=4;
  fprintf(stderr,"%s: max number of points: %5i, using %5i as nx\n",
	  argv[0],nmax,nx);
  // allocate interpolation arrays
  xi=(double *)malloc(nx*sizeof(double));
  yi=(double *)malloc(nx*sizeof(double));
  yis=(double *)malloc(nx*sizeof(double));
  if(!xi || !yi || !yis){
    fprintf(stderr,"%s: memory error: nx: %i\n",argv[0],nx);
    exit(-1);
  }
  /*

    determine extrema of x and the spacing dx based on nx

  */
  xmin=1e20;xmax=-1e20;
  for(i=0;i<nrp;i++){
    if(x[i]<xmin)xmin=x[i];
    if(x[i]>xmax)xmax=x[i];
  }
  dx=(xmax-xmin)/(double)(nx-1);
  fprintf(stderr,"%s: determined xmin: %g xmax: %g, resulting in dx: %g at nx: %i\n",
	  argv[0],xmin,xmax,dx,nx);
  //
  // obtain x values of interpolation 
  // and set the yi (for the mans) and the yis (for the std) to 
  // zero
  //
  for(tx=xmin,i=0;i<nx;i++,tx+=dx){
    xi[i] = tx;
    yi[i] = yis[i] = 0.0;
  }
  /*
    
    loop through x-y pair vectors based on the files to 
    interpolate on the xi[]

  */
  for(nrp=0,i=0;i<nfile;i++){
#if (AVG_ID_INTERPOLATION == 0 )// splines
    // get interpolation weight array
    y2=(double *)malloc(sizeof(double)*n[i]);
    if(!y2){
      fprintf(stderr,"%s: memory error: n[i]: %i\n",argv[0],n[i]);
      exit(-1);
    }
    // obtain spline weights for this file
    spline((x+nrp-1),(y+nrp-1),n[i],1e30,1e30,(y2-1));
    // obtain values of interpolation at the xi
    for(j=0;j<nx;j++){
      splint((x+nrp-1),(y+nrp-1),(y2-1),n[i],xi[j],&ytmp);
      yi[j] += ytmp;
    }
    free(y2);
#else
    for(j=0;j<nx;j++){
      /*
	
	use a polynomial interpolation of order m = AVG_ID_INTERPOLATION

      */
      // find xi in x vector, index will be 1..n style
      // from nrp ... nrp+n[i]
      locate((x+nrp-1),n[i],xi[j],&index);
      // this will still be 1 .. n style
      index = nrp+IMIN(IMAX(index - (order-1)/2,1),n[i] - order + 1);
      // do polynomial interpolation of order "order"
      polint((x+index-2),(y+index-2),order,xi[j],&ytmp,&dy);
      yi[j] += ytmp;
    }
#endif
    nrp += n[i];
  }
  /*
    calculate mean values
  */
  if(nfile)
    for(i=0;i<nx;i++)
      yi[i]/=(double)nfile;
  /*

    calculate the standard deviation at each xi,yi

  */
  if(nfile>1){
    for(nrp=0,i=0;i<nfile;i++){
#if (AVG_ID_INTERPOLATION == 0)// splines
      y2=(double *)malloc(sizeof(double)*n[i]);
      if(!y2){
	fprintf(stderr,"%s: memory error: n[i]: %i\n",argv[0],n[i]);
	exit(-1);
      }
      spline((x+nrp-1),(y+nrp-1),n[i],1e30,1e30,(y2-1));
      for(j=0;j<nx;j++){
	splint((x+nrp-1),(y+nrp-1),(y2-1),n[i],xi[j],&ytmp);
	ytmp -= yi[j];
	yis[j] += ytmp * ytmp;// add squared deviation from mean
      }
      free(y2);
#elif (AVG_ID_INTERPOLATION >= 1) // polynomial
      for(j=0;j<nx;j++){
	locate((x+nrp-1),n[i],xi[j],&index);
	index = nrp + IMIN(IMAX(index - (order-1)/2,1),n[i] - order+1);
	polint((x+index-2),(y+index-2),order,xi[j],&ytmp,&dy);
	ytmp -= yi[j];
	yis[j] += ytmp * ytmp;// add squared deviation from mean
      }
#endif
      nrp += n[i];
    }
    j = nfile-1;
    for(i=0;i<nx;i++)
      yis[i] = sqrt(yis[i]/(double)j);
  }
  //
  // output of average values and standard deviation
  //
  if(nfile)
    for(i=0;i<nx;i++){
      fprintf(stdout,"%12g %12g %12g\n",xi[i],yi[i],yis[i]);
    }
  return 0;
}

/*


  numerical recipes routines from here on,
  precision changed to double


 */


void spline(double *x,double *y,int n,double yp1,double ypn,
	    double *y2)
{
  int i,k;
  double p,qn,sig,un,*u;
  u=(double *)malloc((n+2)*sizeof(double));
  if(!u){fprintf(stderr,"memerror in spline\n");exit(-1);}

  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];

  free(u);
}

void splint(double *xa,double *ya,double *y2a,int n,double x,double *y)
{
  int klo,khi,k;
  double h,b,a;

  klo=1;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) {fprintf(stderr,"Bad xa input to routine splint\n");exit(-1);}
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void polint(double *xa,double *ya,int n,double x,double *y,
	    double *dy)
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c,*d;

  dif=fabs(x-xa[1]);
  c=(double *)malloc((n+2)*sizeof(double));
  d=(double *)malloc((n+2)*sizeof(double));

  if(!c || !d){
    fprintf(stderr,"memerror in polint\n");
    exit(-1);
  }
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0){
	fprintf(stderr,"Error in routine polint\n");
	exit(-1);
      }
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free(c);free(d);
}

void locate(double *xx,int n,double x,int *j)
{
  int ju,jm,jl;
  int ascnd;
  
  jl=0;
  ju=n+1;
  ascnd=(xx[n] > xx[1]);
  while (ju-jl > 1) {
    jm=(ju+jl) >> 1;
    if (x > xx[jm] == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  *j=jl;
}

