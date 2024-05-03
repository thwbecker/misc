#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// for numrec
#define NR_END 1
#define FREE_ARG char*
double *vector(int ,int );
void free_vector(double *,int ,int );
void realft(double *,int ,int );
void four1(double *,int ,int );
void correl(double *,double *,int ,double *);
void twofft(double *,double *,double *,double *,int );
void compute_correl(double **, double **, double **, int, int *);

#define ME {fprintf(stderr,"memory error\n");exit(-1);}
#define INVERSE_LN_TWO  1.44269504088896341
#define log2(x)         (log(x) * INVERSE_LN_TWO)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
