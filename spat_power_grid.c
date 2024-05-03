#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define TWOPID 6.2831853071795865
void avevar(float *,unsigned long ,float *,float *);
void period(float [], float [], int, double, double, float [], float [], int, int *, int *, float *);
double *dvector(long ,long );
void free_dvector(double *,long ,long );
#define FREE_ARG char*
#define NR_END 1
double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define nrerror(x) {fprintf(stderr,"%s\n",x);exit(-1);}
/* 

compute periodograms

*/
int main(int argc, char **argv)
{
  int nx,ny,np,jmax,nout,i,j;
  float *x,y,*z,*px,*py,osample=8.0,prob,hifac=0.25;
  if(argc < 2){
    fprintf(stderr,"%s nx [osample, 5] [hifac, 0.25]\n",argv[0]);
    exit(-1);
  }
  sscanf(argv[1],"%i",&nx);
  if(argc>2)
    sscanf(argv[2],"%f",&osample);
  if(argc>3)
    sscanf(argv[3],"%f",&hifac);
  /*  */
  x=(float *)malloc(sizeof(float)*nx);
  z=(float *)malloc(sizeof(float)*nx);
  /*  */
  np = nx/osample*100;
  px=(float *)malloc(sizeof(float)*np);
  py=(float *)malloc(sizeof(float)*np);
  if(!x || !z || !px || !py){fprintf(stderr,"memerro\n");exit(-1);}

  fprintf(stderr,"%s: expecting data triples with %i entries per row, osample: %g hifac: %g\n",
	  argv[0],nx,osample,hifac);

  i=ny=0;
  while(fscanf(stdin,"%f %f %f",(x+i),&y,(z+i))==3){
    i++;
    if(i == nx){		/* compute periodgram for this row */
      period((x-1),(z-1),nx,osample,hifac,(px-1),(py-1),np,&nout,&jmax,&prob);
      for(j=0;j<nout;j++)
	fprintf(stdout,"%g %g %g\n",px[j],y,log10(py[j]));
      //fprintf(stderr,"%s: compute periodogram %i at %g\n",argv[0],ny,y);
      //fprintf(stdout,"\n\n");
      ny++;
      i=0;
    }
  }
  fprintf(stderr,"%s: computed %i periodograms\n",argv[0],ny);
  return 0;
}


void period(x,y,n,ofac,hifac,px,py,np,nout,jmax,prob)
float *prob,hifac,ofac,px[],py[],x[],y[];
int *jmax,*nout,n,np;
{
	int i,j;
	float ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh,
		sumsy,swtau,var,wtau,xave,xdif,xmax,xmin,yy;
	double arg,wtemp,*wi,*wpi,*wpr,*wr;

	wi=dvector(1,n);
	wpi=dvector(1,n);
	wpr=dvector(1,n);
	wr=dvector(1,n);
	*nout=0.5*ofac*hifac*n;
	if (*nout > np) {
	  fprintf(stderr,"output arrays too short in period, nout: %i np: %i\n",
		  *nout,np);
	  exit(-1);
	}
	avevar(y,n,&ave,&var);
	xmax=xmin=x[1];
	for (j=1;j<=n;j++) {
		if (x[j] > xmax) xmax=x[j];
		if (x[j] < xmin) xmin=x[j];
	}
	xdif=xmax-xmin;
	xave=0.5*(xmax+xmin);
	pymax=0.0;
	pnow=1.0/(xdif*ofac);
	for (j=1;j<=n;j++) {
		arg=TWOPID*((x[j]-xave)*pnow);
		wpr[j] = -2.0*SQR(sin(0.5*arg));
		wpi[j]=sin(arg);
		wr[j]=cos(arg);
		wi[j]=wpi[j];
	}
	for (i=1;i<=(*nout);i++) {
		px[i]=pnow;
		sumsh=sumc=0.0;
		for (j=1;j<=n;j++) {
			c=wr[j];
			s=wi[j];
			sumsh += s*c;
			sumc += (c-s)*(c+s);
		}
		wtau=0.5*atan2(2.0*sumsh,sumc);
		swtau=sin(wtau);
		cwtau=cos(wtau);
		sums=sumc=sumsy=sumcy=0.0;
		for (j=1;j<=n;j++) {
			s=wi[j];
			c=wr[j];
			ss=s*cwtau-c*swtau;
			cc=c*cwtau+s*swtau;
			sums += ss*ss;
			sumc += cc*cc;
			yy=y[j]-ave;
			sumsy += yy*ss;
			sumcy += yy*cc;
			wr[j]=((wtemp=wr[j])*wpr[j]-wi[j]*wpi[j])+wr[j];
			wi[j]=(wi[j]*wpr[j]+wtemp*wpi[j])+wi[j];
		}
		py[i]=0.5*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/var;
		if (py[i] >= pymax) pymax=py[(*jmax=i)];
		pnow += 1.0/(ofac*xdif);
	}
	expy=exp(-pymax);
	effm=2.0*(*nout)/ofac;
	*prob=effm*expy;
	if (*prob > 0.01) *prob=1.0-pow(1.0-expy,effm);
	free_dvector(wr,1,n);
	free_dvector(wpr,1,n);
	free_dvector(wpi,1,n);
	free_dvector(wi,1,n);
}
#undef TWOPID
void avevar(float *data,unsigned long n,float *ave, float *var)
{
	unsigned long j;
	float s,ep;

	for (*ave=0.0,j=1;j<=n;j++) *ave += data[j];
	*ave /= n;
	*var=ep=0.0;
	for (j=1;j<=n;j++) {
		s=data[j]-(*ave);
		ep += s;
		*var += s*s;
	}
	*var=(*var-ep*ep/n)/(n-1);
}

void free_dvector(v,nl,nh)
double *v;
long nh,nl;
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
double *dvector(nl,nh)
long nh,nl;
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}
