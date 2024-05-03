#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/* gnu */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>


struct data {
  int n;
  double * xd;
  double * y;
  double * sigma;
};

int lm_solver(int, int, struct data *d, double *, int);
void print_state(int, gsl_multifit_fdfsolver *);
int erf_f(const gsl_vector *, void *, gsl_vector *);
int erf_df(const gsl_vector *, void *, gsl_matrix *);
int erf_fdf(const gsl_vector *, void *, gsl_vector *, gsl_matrix *);


#define N 40

int main (int argc, char **argv)
{
  unsigned int i;
  const int n = N;
  const int p = 3;
  FILE *out;
  double y[N], sigma[N],xd[N];
  struct data d = { n, xd, y, sigma};
  double solution[3] = { 0, -5.0, 5.0 };
  const gsl_rng_type * type;
  gsl_rng * r;
  gsl_rng_env_setup();
  type = gsl_rng_default;
  r = gsl_rng_alloc (type);
  
  /* This is the data to be fitted */
  out = fopen("tmp.dat","w");
  for (i = 0; i < n; i++)
    {
      xd[i] = i + gsl_ran_gaussian(r,0.1);
      /* d = -1, b = 5, c = 10 */
      y[i] = -1 - 5 * erf(xd[i]/10)  + gsl_ran_gaussian (r, 0.01);
      sigma[i] = .1;
      fprintf (out,"%g %g %g\n", xd[i], y[i], sigma[i]);
    };
  fclose(out);
  gsl_rng_free (r);

  lm_solver(n,p,&d,solution,1);
  fprintf(stderr,"%s: solution: %g %g %g\n",argv[0],solution[0],solution[1],solution[2]);
  return 0;
}

int lm_solver(int n,int p,  struct data *d, double *sol, int verbose)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  gsl_multifit_function_fdf f;
  gsl_vector_view x = gsl_vector_view_array (sol, p);
  int i, iter=0;

  
  f.f = &erf_f;f.df = &erf_df;f.fdf = &erf_fdf; /* functions */
  f.n = n;f.p = p;f.params = d;


  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  if(verbose) 
    print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);
      if(verbose){
	printf ("status = %s\n", gsl_strerror (status));
	
	print_state (iter, s);
      }

      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x,1e-9, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 500);

  gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  { 
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
    if(verbose){
      printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
      
      printf ("d      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
      printf ("b      = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
      printf ("c      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
    }
  }
  if(verbose)
    printf ("status = %s\n", gsl_strerror (status));

  for(i=0;i<p;i++)		/* assign solution */
    sol[i] = gsl_vector_get(s->x, i);

  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  
  return 0;
}

void
print_state (int iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3i x = % 15.8f % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2), 
          gsl_blas_dnrm2 (s->f));
}

int erf_f (const gsl_vector *x, void *data, gsl_vector *f)
{
  
  int n = ((struct data *)data)->n;
  double *xd = ((struct data *)data)->xd;
  double *y = ((struct data *)data)->y;
  double *sigma = ((struct data *) data)->sigma;

  double d = gsl_vector_get (x, 0);
  double b = gsl_vector_get (x, 1);
  double c = gsl_vector_get (x, 2);
  double Yi;
  
  int i;

  for (i = 0; i < n; i++)
    {
      /* Model Yi = d - b erf(x/c) */
      Yi = d - b * erf(xd[i]/c);
      gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
    }

  return GSL_SUCCESS;
}

int erf_df (const gsl_vector *par, void *data, gsl_matrix * J)
{
  const double pi = acos(-1.0);
  const double sqrt_pi = sqrt(pi);
  
  int n = ((struct data *)data)->n;
  double *sigma = ((struct data *) data)->sigma;
  double *xd = ((struct data *) data)->xd;

  double b = gsl_vector_get (par, 1);
  double c = gsl_vector_get (par, 2);

  double fac = -(2*b)/(c*sqrt_pi);

  int i;

  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = d -b * erf(x/c)  */
      /* and the xj are the parameters (d,b,c) */
      double s = sigma[i];
      double xoc = xd[i]/c;
      gsl_matrix_set (J, i, 0,                1 /s); /* df/dd / sigma*/
      gsl_matrix_set (J, i, 1,    -erf(xd[i]/c) /s); /* df/db / sigma*/
      gsl_matrix_set (J, i, 2,fac*exp(-xoc*xoc) /s); /* df/dc / sigma*/

    }
  return GSL_SUCCESS;
}

int
erf_fdf (const gsl_vector * x, void *data,
	 gsl_vector * f, gsl_matrix * J)
{
  erf_f (x, data, f);
  erf_df (x, data, J);

  return GSL_SUCCESS;
}

