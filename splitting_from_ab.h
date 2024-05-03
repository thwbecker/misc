#include "hc.h"
#include "ggrd_grdtrack_util.h"

struct mod{
  double top,bottom;
  struct ggrd_gt a[2];
  int nz;
  struct prem_model prem[1];
};


void read_ab(struct mod *, char *);
void compute_sks(struct mod *, double, double,
		 int , double *, double *,
		 double *, double *,double *,hc_boolean);
