/* pca_util.c */
double calc_mean(double *, int);
void a_equals_b_vector(double *, double *, int);
void add_const_to_vector(double *, double, int);
void compute_cov(double *, double *, int, int);
void compute_aat(double *, double *, int, int);
double vecdotp(double *, double *, int);
FILE *myopen(char *, char *, char *);
