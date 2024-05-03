/* mglin.c */
void mglin(double **, int, int);
void rstrct(double **, double **, int);
void interp(double **, double **, int);
void addint(double **, double **, double **, int);
void slvsml(double **, double **);
void relax(double **, double **, int);
void resid(double **, double **, double **, int);
void copy(double **, double **, int);
void fill0(double **, int);

void sor(double **,double **,double **,double **,double **,double **,
	 double **,int ,double );
void sor_hom(double **,double **,double ,double ,double ,double ,
	     double ,int,double);

