#ifndef COMP_PRECISION
#define COMP_PRECISION double
#endif


#ifndef COMP_EPS
#define COMP_EPS 5.0e-15
#endif

void get_z_rot(COMP_PRECISION [3][3],COMP_PRECISION );
void get_gen_rot(COMP_PRECISION [3][3],COMP_PRECISION,
		 COMP_PRECISION,COMP_PRECISION );
void rotate_vec_2d_angle(COMP_PRECISION *, COMP_PRECISION *,
			 COMP_PRECISION , COMP_PRECISION );
void rotate_vec_3d(COMP_PRECISION [3],COMP_PRECISION [3],
		   COMP_PRECISION [3][3]);

void rotate_mat(COMP_PRECISION [3][3],COMP_PRECISION [3][3],
		COMP_PRECISION [3][3]);

void my_sincos(COMP_PRECISION , COMP_PRECISION *, COMP_PRECISION *);

void calc_polar_basis_at_r(COMP_PRECISION [3][3],COMP_PRECISION *);
void polar_to_cart_mat_at_r3x3(COMP_PRECISION [3][3],COMP_PRECISION [3][3],
			       COMP_PRECISION *);
