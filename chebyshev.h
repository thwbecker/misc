COMP_PRECISION chebev(COMP_PRECISION ,COMP_PRECISION ,
		      COMP_PRECISION *,int ,COMP_PRECISION );
COMP_PRECISION norm_chebev(COMP_PRECISION ,COMP_PRECISION ,
			   COMP_PRECISION *,int ,COMP_PRECISION );
COMP_PRECISION chebnorm(int);
COMP_PRECISION cheb_norm_damping(int,int ,int );
void chder(COMP_PRECISION ,COMP_PRECISION ,
	   COMP_PRECISION *,COMP_PRECISION *,
	   int );
COMP_PRECISION chebev_der_int(COMP_PRECISION,COMP_PRECISION,int);
COMP_PRECISION norm_chebev_der_int(COMP_PRECISION 
				   ,COMP_PRECISION ,int );
