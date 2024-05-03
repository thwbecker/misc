#define COMP_PRECISION double
#include "mglin.h"
#include "nr_defines.h"
#define NPRE 1
#define NPOST 1
#define NGMAX 15

void mglin(double **u, int n,int ncycle)
{
  unsigned int j,jcycle,jj,jpost,jpre,nf,ng=0,ngrid,nn;
  double **ires[NGMAX+1],**irho[NGMAX+1],**irhs[NGMAX+1],**iu[NGMAX+1];
  
  nn=n;
  while (nn >>= 1) ng++;
  if (n != 1+(1L << ng)) nr_error("n-1 must be a power of 2 in mglin.");
  if (ng > NGMAX) nr_error("increase NGMAX in mglin.");
  nn=n/2+1;
  ngrid=ng-1;
  irho[ngrid]=nr_matrix(1,nn,1,nn);
  rstrct(irho[ngrid],u,nn);
  while (nn > 3) {
    nn=nn/2+1;
    irho[--ngrid]=nr_matrix(1,nn,1,nn);
    rstrct(irho[ngrid],irho[ngrid+1],nn);
  }
  nn=3;
  iu[1]=nr_matrix(1,nn,1,nn);
  irhs[1]=nr_matrix(1,nn,1,nn);
  slvsml(iu[1],irho[1]);
  nr_free_matrix(irho[1],1,nn,1,nn);
  ngrid=ng;
  for (j=2;j<=ngrid;j++) {
    nn=2*nn-1;
    iu[j]=nr_matrix(1,nn,1,nn);
    irhs[j]=nr_matrix(1,nn,1,nn);
    ires[j]=nr_matrix(1,nn,1,nn);
    interp(iu[j],iu[j-1],nn);
    copy(irhs[j],(j != ngrid ? irho[j] : u),nn);
    for (jcycle=1;jcycle<=ncycle;jcycle++) {
      nf=nn;
      for (jj=j;jj>=2;jj--) {
	for (jpre=1;jpre<=NPRE;jpre++)
	  relax(iu[jj],irhs[jj],nf);
	resid(ires[jj],iu[jj],irhs[jj],nf);
	nf=nf/2+1;
	rstrct(irhs[jj-1],ires[jj],nf);
	fill0(iu[jj-1],nf);
      }
      slvsml(iu[1],irhs[1]);
      nf=3;
      for (jj=2;jj<=j;jj++) {
	nf=2*nf-1;
	addint(iu[jj],iu[jj-1],ires[jj],nf);
	for (jpost=1;jpost<=NPOST;jpost++)
	  relax(iu[jj],irhs[jj],nf);
      }
    }
  }
  copy(u,iu[ngrid],n);
  for (nn=n,j=ng;j>=2;j--,nn=nn/2+1) {
    nr_free_matrix(ires[j],1,nn,1,nn);
    nr_free_matrix(irhs[j],1,nn,1,nn);
    nr_free_matrix(iu[j],1,nn,1,nn);
    if (j != ng) nr_free_matrix(irho[j],1,nn,1,nn);
  }
  nr_free_matrix(irhs[1],1,3,1,3);
  nr_free_matrix(iu[1],1,3,1,3);
}
#undef NPRE
#undef NPOST
#undef NGMAX

void rstrct(double **uc,double **uf,int nc)
{
  int ic,iif,jc,jf,ncc=2*nc-1;
  
  for (jf=3,jc=2;jc<nc;jc++,jf+=2) {
    for (iif=3,ic=2;ic<nc;ic++,iif+=2) {
      uc[ic][jc]=0.5*uf[iif][jf]+0.125*(uf[iif+1][jf]+uf[iif-1][jf]
					+uf[iif][jf+1]+uf[iif][jf-1]);
    }
  }
  for (jc=1,ic=1;ic<=nc;ic++,jc+=2) {
    uc[ic][1]=uf[jc][1];
    uc[ic][nc]=uf[jc][ncc];
  }
  for (jc=1,ic=1;ic<=nc;ic++,jc+=2) {
    uc[1][ic]=uf[1][jc];
    uc[nc][ic]=uf[ncc][jc];
  }
}

void interp(double **uf,double **uc,int nf)
{
  int ic,iif,jc,jf,nc;
  nc=nf/2+1;
  for (jc=1,jf=1;jc<=nc;jc++,jf+=2)
    for (ic=1;ic<=nc;ic++) uf[2*ic-1][jf]=uc[ic][jc];
  for (jf=1;jf<=nf;jf+=2)
    for (iif=2;iif<nf;iif+=2)
      uf[iif][jf]=0.5*(uf[iif+1][jf]+uf[iif-1][jf]);
  
  for (jf=2;jf<nf;jf+=2)
    for (iif=1;iif <= nf;iif++)
      uf[iif][jf]=0.5*(uf[iif][jf+1]+uf[iif][jf-1]);
}
void addint(double **uf,double **uc,double **res,int nf)
{
  int i,j;
  
  interp(res,uc,nf);
  for (j=1;j<=nf;j++)
    for (i=1;i<=nf;i++)
      uf[i][j] += res[i][j];
}

void slvsml(double **u,double **rhs)
{
  double h=0.5;
  
  fill0(u,3);
  u[2][2] = -h*h*rhs[2][2]/4.0;
}

void relax(double **u,double **rhs,int n)
{
  int i,ipass,isw,j,jsw=1;
  double h,h2;
  
  h=1.0/(n-1);
  h2=h*h;
  for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
    isw=jsw;
    for (j=2;j<n;j++,isw=3-isw)
      for (i=isw+1;i<n;i+=2)
	u[i][j]=0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]
		      +u[i][j-1]-h2*rhs[i][j]);
  }
}
void resid(double **res,double **u,double **rhs,int n)
{
  int i,j;
  double h,h2i;
  
  h=1.0/(n-1);
  h2i=1.0/(h*h);
  for (j=2;j<n;j++)
    for (i=2;i<n;i++)
      res[i][j] = -h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-
			4.0*u[i][j])+rhs[i][j];
  for (i=1;i<=n;i++)
    res[i][1]=res[i][n]=res[1][i]=res[n][i]=0.0;
}

void copy(double **aout,double **ain,int n)
{
  int i,j;
  for (i=1;i<=n;i++)
    for (j=1;j<=n;j++)
      aout[j][i]=ain[j][i];
  
}

void fill0(double **u, int n)
{
  int i,j;
  for (j=1;j<=n;j++)
    for (i=1;i<=n;i++)
      u[i][j]=0.0;
}

/* successive over-relaxation  */
#include <math.h>
#define MAXITS 1000
#define EPS 1.0e-5

void sor(double **a,double **b,double **c,double **d,double **e,double **f,
	 double **u,int jmax,double rjac)

{
  int ipass,j,jsw,l,lsw,n;
  double anorm,anormf=0.0,omega=1.0,resid;
  
  for (j=2;j<jmax;j++)
    for (l=2;l<jmax;l++)
      anormf += fabs(f[j][l]);
  for (n=1;n<=MAXITS;n++) {
    anorm=0.0;
    jsw=1;
    for (ipass=1;ipass<=2;ipass++) {
      lsw=jsw;
      for (j=2;j<jmax;j++) {
	for (l=lsw+1;l<jmax;l+=2) {
	  resid=a[j][l]*u[j+1][l]
	    +b[j][l]*u[j-1][l]
	    +c[j][l]*u[j][l+1]
	    +d[j][l]*u[j][l-1]
	    +e[j][l]*u[j][l]
	    -f[j][l];
	  anorm += fabs(resid);
	  u[j][l] -= omega*resid/e[j][l];
	}
	lsw=3-lsw;
      }
      jsw=3-jsw;
      omega=(n == 1 && ipass == 1 ? 1.0/(1.0-0.5*rjac*rjac) :
	     1.0/(1.0-0.25*rjac*rjac*omega));
    }
    if (anorm < EPS*anormf) return;
  }
  fprintf(stderr,"MAXITS exceeded in SOR\n");
  exit(-1);
}

void sor_hom(double **f,double **u,double a,double b,double c,double d,
	     double e,int jmax,double rjac)

{
  int ipass,j,jsw,l,lsw,n;
  double anorm,anormf=0.0,omega=1.0,resid;
  
  for (j=2;j<jmax;j++)
    for (l=2;l<jmax;l++)
      anormf += fabs(f[j][l]);
  for (n=1;n<=MAXITS;n++) {
    anorm=0.0;
    jsw=1;
    for (ipass=1;ipass<=2;ipass++) {
      lsw=jsw;
      for (j=2;j<jmax;j++) {
	for (l=lsw+1;l<jmax;l+=2) {
	  resid=a*u[j+1][l]
	    +b*u[j-1][l]
	    +c*u[j][l+1]
	    +d*u[j][l-1]
	    +e*u[j][l]
	    -f[j][l];
	  anorm += fabs(resid);
	  u[j][l] -= omega*resid/e;
	}
	lsw=3-lsw;
      }
      jsw=3-jsw;
      omega=(n == 1 && ipass == 1 ? 1.0/(1.0-0.5*rjac*rjac) :
	     1.0/(1.0-0.25*rjac*rjac*omega));
    }
    if (anorm < EPS*anormf) return;
  }
  fprintf(stderr,"MAXITS exceeded in SOR\n");
  exit(-1);


}

#undef MAXITS
#undef EPS


