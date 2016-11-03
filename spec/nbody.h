#ifndef __nbody_h__
#define __nbody_h__
#include <string.h>

/* const real_t G_const = 6.67408e-11; */ /* m^3 kg^{-1} s^{-2} */
/* const real_t G_const = 4.302e-3;*/ /* pc M_sun^{-1}(km/s)^2 */
#define soft_eps 0.1
#define G_const 6.67408e-11
#define theta 0.01
#define box_length 10e12
#define start_t_const 0
#define t_end_const 15778000.0
#define h_const 25000.0
#define input_filename "../data/planets_1.tab" 
/* const int Order = 4;*/ /* for rk4 */
/* const real_t eps = 1e-1;*/ /* this is just a wild guess! */

typedef double real_t; 
typedef struct rvector_s {
    real_t x[3];
} rvector_t;

typedef struct particle_t {
	 real_t pos[3];
	 real_t vel[3];
	 real_t mass;
} particle_t;


typedef struct octree_node_s octree_node_t;

struct octree_node_s{
    real_t length;
    particle_t *particle;
    octree_node_t *paren;
    octree_node_t *children;
};

typedef struct area_s{
    real_t x_min;
    real_t x_max;
    real_t y_min;
    real_t y_max;
    real_t z_min;
    real_t z_max;
} area_t;
    

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
rvector_t nbody_bh_calculate_force(particle_t *particle, octree_node_t *tree);
void nbody_add_particle_particle_interaction(particle_t *particle_i, particle_t *particle_j, rvector_t *res);
void octree_insert_node(particle_t *particle, octree_node_t *tree, rvector_t center,real_t length);
octree_node_t * octree_generate_tree(long N, particle_t *particles, rvector_t center, real_t length);
void octree_free(octree_node_t *tree);
void octree_pretty_print(octree_node_t *tree);
void particle_pretty_print(particle_t *particle);
void RKstep ( rhs_function_t f,  real_t t, real_t *y, real_t *ynew, real_t h, long N, particle_t *particles);
void nbodyprob(real_t t, real_t *y, real_t *x, long N, particle_t *particles);
void print_step(real_t t, real_t *y, char *output_filename_base, long N, particle_t *particles);
#endif
