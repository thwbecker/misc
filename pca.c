#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "precision.h"
#include "nr_defines.h"
#include "misc.h"
#include "pca.h"
extern void rs_(int *, int *, double *,double *, int *, 
		double *, double *,double *, int *);


/* 
   

principal component analysis


reads in columns of data with dim number of rows (dimensions),
performs PCA analysis, and writes eigensystem as well as projected
data to files


$Id: pca.c,v 1.1 2005/07/08 20:38:50 becker Exp becker $


*/


int main(int argc, char **argv)
{
  COMP_PRECISION *d,*dm,*ds,*dtmp,*s,*c,*fv1,*fv2,*ps,*eigenr,*v,comp,*tvec,
    tnorm_s,tdev,dev,*rec,*crms;
  int i,j,k,l,n,dim,ndim,ierr,dimm1,iuse;
  static my_boolean 
    debug_file_io = FALSE, 	/* all kinds of stuff */
    reduced_data_file_io = TRUE, /* data - mean */
    file_io = TRUE,		/* general file I/O? */
    make_positive=FALSE;	/* eigenvectors */
  FILE *out;
  static int rs_matz = 1;
  /* 
     deal with command line arguments 
  */
  if(argc != 2){
    fprintf(stderr,"%s d\nperform principal component analysis of n data with d dimensions\n",
	    argv[0]);
    fprintf(stderr,"\treads in rows with d columns from stdin\n");
    exit(-1);
  }
  sscanf(argv[1],"%i",&dim);
  dimm1 = dim-1;
  fprintf(stderr,"%s: expecting %i dimensional data from stdin\n",
	  argv[0],dim);
  if(debug_file_io){
    reduced_data_file_io = TRUE;
    file_io = TRUE;
  }
  /* data input */
  my_vecalloc(&d,1,argv[0]);
  i=0;
  while(fscanf(stdin,"%lf",(d+i++))==1)
    my_vecrealloc(&d,(i+1),argv[0]);
  if(i%dim == 0){
    fprintf(stderr,"%s: data mismatch, read %i entries\n",argv[0],i);
    exit(-1);
  }
  n = i/dim;			/* number of data points */
  ndim = n * dim;			/* total matrix size */
  /*
    
     resort data so that x is [d][n]

  */
  my_vecalloc(&dtmp,ndim,argv[0]);
  for(i=0;i<dim;i++)
    for(j=0;j<n;j++)
      dtmp[i*n+j] = d[j*dim+i];
  a_equals_b_vector(d,dtmp,ndim);
  free(dtmp);
  fprintf(stderr,"%s: read %i data entries\n",argv[0],n);
  /* 
     compute the means and stddevs
  */
  my_vecalloc(&dm,dim,argv[0]);
  for(i=0;i<dim;i++)
    dm[i] = calc_mean((d+i*n),n);
  /* 
     remove means from data 
  */
  my_vecalloc(&s,ndim,argv[0]);
  for(i=0;i<dim;i++)
    for(j=0;j<n;j++)
      s[i*n+j] = d[i*n+j] - dm[i];
  /* compute total norm */
  for(tnorm_s=0.0,i=0;i<n;i++){
    for(j=0;j < dim;j++)
      tnorm_s += s[j*n+i]*s[j*n+i];
  }
  tnorm_s = sqrt(tnorm_s);
  if(reduced_data_file_io){
    /* 
       write means 
    */
    fprintf(stderr,"%s: writing mean values to pca.mean.dat\n",argv[0]);
    out = myopen("pca.mean.dat","w",argv[0]);
    for(j=0;j < dim;j++)
      fprintf(out,"%lg\t ",dm[j]);
    fprintf(out,"\n");
    fclose(out);
    /* 
       writing s = x - m  data to stdout 
    */
    fprintf(stderr,"%s: writing reduced original data pca.s.dat\n",argv[0]);
    out = myopen("pca.s.dat","w",argv[0]);
    for(i=0;i < n;i++){
      for(j=0;j<dim;j++)
	fprintf(out,"%lg\t",s[j*n+i]);
      fprintf(out,"\n");
    }
    fclose(out);
  }
  /* 
     compute the correlation matrix,
     (s.s^T)/N
  */
  my_vecalloc(&c,dim*dim,argv[0]);
  compute_cov(c,s,dim,n);
  /* 
     compute eigen vectors
  */
  my_vecalloc(&eigenr,dim,argv[0]);my_vecalloc(&fv1,dim,argv[0]);
  my_vecalloc(&fv2,dim,argv[0]); my_vecalloc(&v,dim*dim+10,argv[0]);
  /* 
     returns eigenvalues (eigenr) and normalized vectors sorted from
     smallest to largest (v has the vectors in rows)

  */
  rs_(&dim, &dim, c, eigenr, &rs_matz, v, fv1, fv2, &ierr);
  free(fv1);free(fv2);my_vecrealloc(&v,dim*dim,argv[0]);
  if(ierr){
    fprintf(stderr,"%s: error %i in rs eigen routine\n",
	    argv[0],ierr);
    exit(-1);
  }
  if(make_positive){
    /* make sure the eigenvectors are such that the first dimension is
       positive  */
    for(i=0;i<dim;i++){
      if(v[i*dim+0] < 0)
	for(j=0;j<dim;j++)
	  v[i*dim+j] *= -1.0;
    }
  }
  if(debug_file_io && (dim <= 3)){
    /* 
       output of eigenvectors to gnuplot file 
    */
    fprintf(stderr,"%s: writing gnuplot file to pca.gpl\n",argv[0]);
    out = myopen("pca.gpl","w",argv[0]);
    for(i=0;i<dim;i++){
      if(dim == 2)
	fprintf(out,"set arrow from 0, 0 to ");
      else
	fprintf(out,"set arrow from 0, 0, 0 to ");
      for(comp=0.0,j=0;j<dim;j++){
	comp += v[i*dim+j]*v[i*dim+j];
	if(j==dimm1)
	  fprintf(out,"%11g\n",v[i*dim+j]);
	else
	  fprintf(out,"%11g, ",v[i*dim+j]);
      }
      if(fabs(sqrt(comp) - 1.0)> EPS_COMP_PREC){
	fprintf(stderr,"%s: normalization error: %.15e\n",argv[0],comp);
	exit(-1);
      }
    }
    if(dim == 2)
      fprintf(out,"plot 'pca.s.dat' ps 0.1\n");
    else
      fprintf(out,"splot 'pca.s.dat' ps 0.1\n");
    fclose(out);
  }
  /* 
     project the data into new system
  */
  my_vecalloc(&ps,ndim,argv[0]);
  my_vecalloc(&crms,dim,argv[0]);
  for(i=0;i<dim;i++)crms[i] = 0.0;
  for(i=0;i < n;i++){		/* data loop */
    for(j=0;j<dim;j++){		/* vector loop */
      for(comp=0.0,k=0;k<dim;k++)	/* projected data in new coordinate system */
	comp  += s[k*n+i] * v[j*dim+k];
      ps[j*n+i] = comp;
      crms[j] += comp*comp;
    }
  }
  for(i=0;i<dim;i++)
    crms[i] = sqrt(crms[i]/n);
  
 
  /* test the recovery */
  my_vecalloc(&tvec,dim,argv[0]);
  my_vecalloc(&rec,dim,argv[0]);
  for(i=1;i <= dim;i++){		/* number of factors to use */
    tdev = 0.0;
    for(j=0;j < n;j++){		/* data loop */
      for(k=0;k < dim;k++)	/* clear test vector */
	tvec[k] = 0.0;
      for(iuse=dimm1,k=0;k < i;k++,iuse--){	/* reconstruction loop, iuse
						   is the eigenvector to
						   use */
	for(l=0;l < dim;l++)	/* assemble test vector */
	  tvec[l] += ps[iuse*n+j] * v[iuse*dim+l];
      }	/* end reconstruction */
      /* 
	 compute misfit 
      */
      for(k=0;k < dim;k++){
	comp =  tvec[k] - s[k*n+j];
	tdev += comp*comp;
      }
    } /* end data loop */
    tdev = sqrt(tdev);
    rec[i-1] = 1.0 - tdev/tnorm_s;	/* relative recovery power */
    fprintf(stderr,"%s: relative reconstruction power for %i eigenvalue(s): %.7f\n",
	    argv[0],i,rec[i-1]);
  } /* end recovery test loop */
  if(file_io){
    /* 
       
    general output

    */
    /* 
       output of eigenvalues 
    */
    fprintf(stderr,"%s: writing eigensystem to pca.eigen.dat, format:\n",argv[0]);
    fprintf(stderr,"%s: i eigen_i\t\teigen_i/max_eigen\tinverse_recovery_err ps_component_norm\t vector\n",
	    argv[0]);
    out=myopen("pca.eigen.dat","w",argv[0]); 
    for(i=0;i<dim;i++){
      fprintf(out,"%3i\t%.5e %.7f\t%.7f\t%11g\t",
	      i+1,eigenr[i],eigenr[i]/eigenr[dimm1],
	      rec[dimm1-i],crms[i]);
      for(j=0;j<dim;j++)
	fprintf(out,"%9.6f ",v[i*dim+j]);
      fprintf(out,"\n");
    }

    fclose(out);
    
    fprintf(stderr,"%s: writing projected data to pca.ps.dat\n",argv[0]); 
    out = myopen("pca.ps.dat","w",argv[0]); 
    for(i=0;i<n;i++){ 
      for(j=0;j<dim;j++) 
	fprintf(out,"%11g ",ps[j*n+i]); 
      fprintf(out,"\n"); 
    } 
    fclose(out); 
  } /* end output branch */
  /* free arrays */
  free(eigenr);free(v);free(ps);
  free(d);free(dm);free(tvec);free(rec);
  free(s);free(c);free(crms);
    
}



