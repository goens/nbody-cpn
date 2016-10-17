#include <string.h>

typedef double real_t; 
typedef real_t rvector_t[3];
typedef struct particle_t {
	 rvector_t pos;
	 rvector_t vel;
	 real_t mass;
} particle_t;


/* const real_t G_const = 6.67408e-11; */ /* m^3 kg^{-1} s^{-2} */
const real_t G_const = 4.302e-3; /* pc M_sun^{-1}(km/s)^2

const int Order = 4; //for rk4
const real_t eps = 1e-1; //this is just a wild guess!

/* Improvement: use runge-kutta-fehlberg for an evaluation of order 5.
 * see: http://www.aip.de/groups/soe/local/numres/bookcpdf/c16-2.pdf
 */
typedef void (*rhs_function_t) (real_t t, real_t *y, real_t *x,  long N, particle_t *particles); 
typedef void (*rk_init_function_t) (int *len_Alpha, real_t *Alpha,
												int *len_Beta, int *height_Beta, real_t *Beta,
												int *len_Gamma, real_t *Gamma);

void set_rk4(int *len_Alpha, real_t *Alpha,
				 int *len_Beta, int *height_Beta, real_t *Beta,
				 int *len_Gamma, real_t *Gamma){

	 (*len_Alpha) = 4; 
	 (*len_Gamma) = 4; 
	 (*len_Beta) = 4; 
	 (*height_Beta) = 4; 

	 const real_t Alpha_rk4[] = {0., 1./2., 1./2., 1. };
	 memcpy(Alpha, Alpha_rk4,  sizeof Alpha_rk4);

	 const real_t Beta_rk4[] = { 0., 0., 0., 0., 
								  1./2., 0.,0.,0.,
								  0., 1./2.,0.,0.,
								  0.,0.,1.,0.};
	 memcpy(Beta, Beta_rk4,  sizeof Beta_rk4);

	 const real_t Gamma_rk4[] = {1./6., 1./3., 1./3., 1./6. };
	 memcpy(Gamma, Gamma_rk4,  sizeof Gamma_rk4);

}

void set_euler(int *len_Alpha, real_t *Alpha,
				 int *len_Beta, int *height_Beta, real_t *Beta,
				 int *len_Gamma, real_t *Gamma){

	 (*len_Alpha) = 1; 
	 (*len_Gamma) = 1; 
	 (*len_Beta) = 1; 
	 (*height_Beta) = 1; 

	 const real_t Alpha_euler[] = {0.};
	 memcpy(Alpha, Alpha_euler,  sizeof Alpha_euler);

	 const real_t Beta_euler[] = {0.};
	 memcpy(Beta, Beta_euler,  sizeof Beta_euler);

	 const real_t Gamma_euler[] = {1.};
	 memcpy(Gamma, Gamma_euler,  sizeof Gamma_euler);

}

void RKstep ( rhs_function_t f,  real_t t, real_t *y, real_t *ynew, real_t h, long N, particle_t *particles);
void nbodyprob(real_t t, real_t *y, real_t *x, long N, particle_t *particles);
void print_step(real_t t, real_t *y, char *output_filename_base, long N, particle_t *particles);
