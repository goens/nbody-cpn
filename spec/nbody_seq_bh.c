#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* #include <gsl/gsl_cblas.h> */

#include "nbody.h"

/*
 * DATASET:
 * Dubinksi 1995 : The file dubinski.tab.gz contains a simple (free format) table of 7 columns and 81920 rows:
 * the masses and six phase space coordinates (x,y,z,vx,vy,vz) of 81920 particles that define the initial conditions
 * for a MW/Andromeda collision. 
 * Contact John Dubinski (dubinski@cita.utoronto.ca) for more details 
 * (see also: Dubinski, Mihos & Hernquist, 1996, ApJ, 462, 576). Its the model B collision discussed in this paper). 

 * \delta t = 0.1, this corresponds to 1.8e6 years when scaled to physical units for the galaxy
 * Plummer softening radius \epsilon = 0.025 (80pc) 
 */


/*
 * This is a Runge-Kutta step
 * For a given function f, it takes:
 *   - the time (or in general, ODE variable) t
 *   - the (current) solution vector (approximation) y
 *   - the step length for the next step h
 *
 */

octree_node_t * octree_generate_tree(long N, particle_t *particles, rvector_t center, real_t length){
    octree_node_t *tree = calloc(1,sizeof(octree_node_t));
    long i;
    int l;
    for(i=0;i<N;i++){
        for(l = 0; l <3;l++){
            if( abs(center.x[l] - particles[i].pos[l]) > length){
                printf("Error, particle does not fit into domain (center: %lf %lf %lf, length %lf, particle: ", center.x[0],center.x[1],center.x[2],length);
                particle_pretty_print(particles + i);
                printf(")\n");
                exit(1);
            }
        }
        octree_insert_node(particles + i, tree, center, length);
    }
    return tree;
}

void octree_insert_node(particle_t *particle, octree_node_t *tree, rvector_t center,real_t length){
    int j;
    if(length == 0){
	printf("error, tried to insert a particle into a volumeless cube\n");
	exit(1);
	}


    if(tree->particle == NULL){ /* unocupied octant */
        tree->particle = particle;
    }
    else{
        if(tree->children == NULL){ /* ocupied octant, leaf: repartition */
            tree->children = calloc(8,sizeof(octree_node_t));
            for(j=0;j<8;j++){
                tree->children[j].paren = tree;
                tree->children[j].length = length; /* full length of octant is double the distance from center to edge */
            }

            /*        --------------                 */
            /*      / |  3  |  7  /|                 */
            /*     /  -----------/--                 */
            /*    /   |  2  |  6/  |                 */
            /*   /   ----------/---   z              */
            /*   --------------   /    ^  ^y         */
            /*   |  1  |  5   |  /     | /           */
            /*   -------------- /      |/            */
            /*   |  0  |  4   |/        --> x        */
            /*   --------------                      */

            /* decide which octant:                  */
            real_t new_length = length / 2.;
            /* printf("length = %lf; new length = %lf\n", length, new_length); */
            int octant_index = 4*(particle->pos[0] < center.x[0]) + 2*(particle->pos[1] < center.x[1]) + (particle->pos[2] < center.x[2]);
            rvector_t new_center;
            for(j=0;j<3;j++){
                /* printf(" {delta = (%lf - %lf)*%lf}", (real_t)(particle->pos[j] < center.x[j]), 1/2. , new_length); */
                new_center.x[j] = center.x[j] + ((real_t)(particle->pos[j] > center.x[j]) - 1/2.) * new_length;
                /* printf(" %lf ",new_center.x[j]); */
            }
            /* printf("\n"); */

            octree_insert_node(particle,&(tree->children[octant_index]),new_center,new_length);

            /* /\* check that particles are not the exactly at the same position *\/ */
            /* /\* this is going to be a problem with cpn, better move it to the nbody step *\/ */
            /* int equal = 1; */
            /* for(j=0;j<3;j++){ */
            /*     equal = equal && (tree->particle->pos[j] == particle->pos[j]); */
            /* } */
            /* if(equal){ */
            /*     for(j=0;j<3;j++){ */
            /*         tree->particle->pos[j] += soft_eps; */
            /*         particle->pos[j] -= soft_eps; */
            /*     } */
            /* } */
                

            /* also for the existing particle */
            octant_index = 4*(tree->particle->pos[0] < center.x[0]) + 2*(tree->particle->pos[1] < center.x[1]) + (tree->particle->pos[2] < center.x[2]);
            /* printf("new center for p = ["); */
            /* particle_pretty_print(tree->particle); */
            /* printf("] = "); */
            for(j=0;j<3;j++){
                new_center.x[j] = center.x[j] + ((real_t)(tree->particle->pos[j] > center.x[j]) - 1/2.) * new_length;
                /* printf(" %lf ",new_center.x[j]); */
            }
            /* printf("\n"); */

            octree_insert_node(tree->particle,&(tree->children[octant_index]),new_center,new_length);

            /* create a center of mass particle */
            particle_t * center_of_mass = calloc(1,sizeof(particle_t));
            center_of_mass->mass = particle->mass + tree->particle->mass;
            for(j=0;j<3;j++){
                center_of_mass->pos[j] = (particle-> mass * particle->pos[j] + tree->particle->mass * tree->particle->pos[j])/center_of_mass->mass;
                center_of_mass->vel[j] = 0; /* init, but this value has no meaning */
            }
            tree->particle = center_of_mass;
        }
        else{ /* tree has children, just update center of mass */
            for(j=0;j<3;j++){
                tree->particle->pos[j] = (tree->particle->mass * tree->particle->pos[j] + particle->mass*particle->pos[j])/ (tree->particle->mass + particle->mass) ;
            }
            tree->particle->mass += particle->mass;

            /* insert in subtree */
            real_t new_length = length / 2.;
            int octant_index = 4*(particle->pos[0] < center.x[0]) + 2*(particle->pos[1] < center.x[1]) + (particle->pos[2] < center.x[2]);
            rvector_t new_center;
            for(j=0;j<3;j++){
                new_center.x[j] = center.x[j] + ((particle->pos[j] > center.x[j]) - 1/2.) * new_length;
            }
            octree_insert_node(particle,&(tree->children[octant_index]),new_center,new_length);
        }

    }
}

void octree_free(octree_node_t *tree){
    int j;
    if(tree == NULL) return;

    if(tree->particle != NULL && tree->children != NULL){ 
        /* particle is a center of mass particle, free it */
        free(tree->particle); 
        tree->particle = NULL;
    }

    if(tree->children != NULL){
        for(j=0;j<8;j++){
            octree_free(tree->children+j);
        }
        free(tree->children);
        tree->children=NULL;
    }

    /* free root node */
    if(tree->paren == NULL) free(tree);
}

void octree_pretty_print(octree_node_t *tree){
    int j = 0;

    if( tree == NULL ) return;
    printf(" (");
    if( tree->particle != NULL){
        printf(" [");
        particle_pretty_print(tree->particle);
        printf("] ");
    }
    if(tree->children != NULL){
        for(j=0;j<8;j++){
            octree_pretty_print(&(tree->children[j]));
        }
    }
    printf(")");

}

void particle_pretty_print(particle_t *particle){
    printf("%lf %lf %lf %lf %lf %lf %lf", particle->mass, particle->pos[0], particle->pos[1], particle->pos[2],particle->vel[0],particle->vel[1],particle->vel[2]);
}

int read_particle(FILE *fp, particle_t *res){
	 int ret = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf\n", &(res->mass), &(res->pos[0]), &(res->pos[1]), &(res->pos[2]), &(res->vel[0]), &(res->vel[0]), &(res->vel[0]));  
    /* testing: no initial velocities */
	 res->vel[0] = 0;
	 res->vel[1] = 0;
	 res->vel[2] = 0;
	 
	 return ret;
}


void read_file(FILE* fp, particle_t *particles){
	 long i = 0;
	 while(read_particle(fp, &(particles[i]) ) == 7) i++;
}


void RKstep( rhs_function_t f, real_t t, particle_t *particles_out, particle_t *particles_in, real_t h, long N){
  
  long i,j,l;
  /* indices for: 
     i: equations (particles)
     l: equations (components)
	  j: RK Butcher-Tableau
  */

  /*initialize with the maximal lengths from the different methods */
  int m = 4;
  const int n = 6; //components (2*dimensions)
  
  /* Butcher Tableau */
  int len_Alpha, len_Gamma, len_Beta, height_Beta;
  real_t Alpha[m];
  real_t Gamma[m];
  real_t Beta[m*m];
  set_rk4(&len_Alpha, Alpha, &len_Beta, &height_Beta, Beta, &len_Gamma, Gamma);
  m = len_Alpha;
  long k;

  /* adaptive step parameters */
  /* rvector_t *e; */
  /* real_t eps2,c; */

  particle_t *part_new = calloc(N, sizeof(particle_t));
  particle_t *sum =  calloc(N, sizeof(particle_t));
  particle_t **K = calloc(m,sizeof(particle_t *));
  for(j=0;j<m;j++){
		K[j] = calloc(N,sizeof(particle_t));
  }/* like: K[m][n*N];  this way K[0] ... K[m-1] are arrays of length N */

  /* do{ */ /*this one is for the adaptive version */
	 
	 /* execute the ODE function m times */
	 for(j=0;j<m;j++){

		  for(i=0;i<N;i++)
				for(l=0;l<n/2;l++){
					 part_new[i].pos[l] = particles_in[i].pos[l]; 
					 part_new[i].vel[l] = particles_in[i].vel[l]; 
		  }
		  /* update the arguments to calculate the next K[j] */
		  for(i=0;i<N;i++){

		  
				for(l=0;l<3;l++){
					 sum[i].pos[l] = 0;
					 sum[i].vel[l] = 0;

					 for(k=0;k<j;k++){
						  sum[i].pos[l] += Beta[j*m+k]*(K[k]+i)->pos[l];
						  sum[i].vel[l] += Beta[j*m+k]*K[k][i].vel[l];
					 }
					 part_new[i].pos[l] += h * sum[i].pos[l];
					 part_new[i].vel[l] += h * sum[i].vel[l];
						  
				}
		  }
		  (*f)(t + Alpha[j] * h, part_new,(particle_t *) &K[j][0], N); /* K[j] = (*f)(t + Alpha[j] * h, ynew + (h * sum)); */
		  /*printf("finished. K[%ld][0] = %lf\n",j,K[j][0]); */
	 } 

 
	 for(i=0;i<N;i++){

	   for(l=0;l<3;l++){
			 sum[i].pos[l] = 0;
			 sum[i].vel[l] = 0;

		  for(j=0;j<m;j++){
         /* if(i==0) */
         /*     printf("sum[%ld] = %lf*%lf\n",l,Gamma[j],K[j][l]); */
				sum[i].pos[l]+=Gamma[j]*K[j][i].pos[l]; //could rewrite in atlas
		  }

		   part_new[i].pos[l] = particles_in[i].pos[l] + h*sum[i].pos[l]; //atlas?
		   part_new[i].vel[l] = particles_in[i].vel[l] + h*sum[i].vel[l]; //atlas?
         /* if(i==0) */
         /*     printf("ynew[%ld] = %lf + %lf*%lf\n",l,y[l],h,sum[l]); */
	   }

	 }
	 /* error estimation (for the adaptive version) */

	 /* for(j=0;j<n;j++) sum[j] = 0; */

	 /* for(j=0;j<m;j++){ */
	 /* 	  sum+=Delta(j)*K(j); //atlas */
	 /* } */
	 /* e = h*sum; */
	 /* eps2 = e.NormMax(); //atlas? */
	 /* //cout << "eps = " << eps << " eps2 = " << eps2;   */
	 /* c = (real)0.9 * powl(eps/eps2  ,(1./(1.+Order))); */
	 /* t = t+h; */

	 /* if(eps2>eps){ */
	 /* 	t -= h; */
	 /* } */

	 /* if(c>5) c = 5; */
	 /* if(c<0.1) c = 0.1; */
	 /* h *= c; */
	 
	 
	 /* }while(eps2>eps); */


	 /* update results after calculating */
	 for(i=0;i<N;i++)
		  for(l=0;l<n/2;l++){
				particles_out[i].pos[l] = part_new[i].pos[l];
				particles_out[i].vel[l] = part_new[i].vel[l];
		  }
	 
	 /* clean up */

	 for(j=0;j<m;j++){
		  free(K[j]);
	 }
	 free(part_new);
	 free(sum);
	 free(K);

}

/* explicitly not inplace to be able to convert it easier to cpn */
void nbodyprob(real_t t, particle_t *particles_in, particle_t *particles_out, long N){
    long i,l;

    const int n = 6; /*components (2*dimensions) */
    /* printf("executing nbody function with t = %lf",t); */
    rvector_t zero = {{ 0., 0., 0.}};
    rvector_t f_i;
    octree_node_t *tree = octree_generate_tree(N, particles_in, zero,box_length);
    /*octree_pretty_print(tree); */
    /*printf("\n"); */
    
  
    for(i=0;i<N;i++){
        for(l=0;l<n/2;l++) /* the first 3 components */
            particles_out[i].pos[l] = particles_in[i].vel[l];  
        for(l=n*i+n/2;l<n*(i+1);l++) /* the second 3 components (where f will be applied) */
            particles_out[i].vel[l] = 0;
        /* if (i==0) printf( "x = %lf. ", x[0]); */
        f_i = nbody_bh_calculate_force(particles_in+i,tree);
		  /* printf("f[%ld]=(%lf,%lf,%lf)",i,f_i.x[0],f_i.x[1],f_i.x[2]); */
        for(l=0;l<n/2;l++){
            particles_out[i].vel[l] = f_i.x[l];
        }
    }

    /* octree_free(tree); */
        /* if (i== 0) printf( " (af. int.) x[0] = %lf. ", x[0]); */
}

void nbody_add_particle_particle_interaction(particle_t *particle_i, particle_t *particle_j, rvector_t *res){
    real_t norm;
    int l;
	 real_t soft_eps_sq = soft_eps * soft_eps; 
    const int n = 6;
    norm = 0;
    for(l = 0; l<n/2;l++){
        norm += pow(particle_j->pos[l]-particle_i->pos[l],2);
        if(norm < soft_eps_sq) /* minimal distance to avoid numeric errors! */
            norm = soft_eps_sq;
    }
    norm = pow(norm,3./2.);
    for(l =0; l<n/2; l++){
        res->x[l] = (particle_j->mass * G_const/norm * (particle_j->pos[l]-particle_i->pos[l]));
    }
}

rvector_t nbody_bh_calculate_force(particle_t *particle, octree_node_t *tree){
    rvector_t f = {{0., 0., 0.}};
    int i,l;
    const int n = 6;

    if(tree->particle == NULL)
        return f;

    if(tree->children == NULL){ /* external node */
        if(tree->particle == particle){
            return f; /*  no interaction with itself */
        }
        else{ /* calculate particle-particle interaction */
            nbody_add_particle_particle_interaction(tree->particle,particle,&f);
        }
    }
    else{ /* not external node. see it it should be considered more */
        real_t ratio = tree->length;
        real_t diff = 0;
        for(l=0;l<n/2;l++)
            diff += pow( tree->particle->pos[l] - particle->pos[l], 2);
        ratio = ratio / sqrt(diff);

        if(ratio < theta){ /* sufficiently far away, approximate with center of mass */
            nbody_add_particle_particle_interaction(tree->particle,particle,&f);
        }
        else{ /* not far away enough, add the component from all children */
            for(i=0;i<8;i++){
                rvector_t temp;
                temp = nbody_bh_calculate_force(particle,tree->children + i);
                for(l=0;l<3;l++){
                    f.x[l] += temp.x[l];

                }
            }
        }

    }
    return f;
}
void print_step(real_t t, char *filename_base, long N, particle_t *particles){
  long i;
  char output_string[strlen(filename_base) + 15];
  strcpy(output_string,filename_base);
  char time[15];
  sprintf(time,".csv.%010.0lf",t*1e6);
  strcat(output_string,time);
  struct stat st = {0};

  if (stat("output", &st) == -1) {
		mkdir("output", 0700);
  }
  FILE *output_file = fopen(output_string,"w");

  /* fprintf(output_file,"time t = %lf: \n", t); */
	 for(i=0;i<N;i++){
        fprintf(output_file,"%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",particles[i].mass,particles[i].pos[0],particles[i].pos[1],particles[i].pos[2],particles[i].vel[0],particles[i].vel[1],particles[i].vel[2]);
        /* fwrite(y + (6*i),sizeof(real_t),3,output_file); */
        /* fwrite(&particles[i].mass,sizeof(real_t),1,output_file); */
	 }
  fclose(output_file);
}

int main(){


	 /****************************/
	 /* read particles from file */
	 /****************************/

	 char *name = input_filename;
	 char *output_file = "output/nbody_sim";
	 FILE *fp = fopen(name,"r");
	 
	 long N = 0;
	 if(fp == NULL){
		  printf("An error ocurred while reading the file: %s\n", name);
	 }
	 /* a bit inefficient: determine the number of lines first for allocation */
	 while(!feof(fp))
	 {
		  char ch = fgetc(fp);
		  if(ch == '\n')
		  {
				N++;
		  }
	 }
	 rewind(fp);

	 particle_t *particles = calloc(N,sizeof(particle_t));
	 read_file(fp,particles);
	 printf("read a file with %ld particles\n",N);

	 /*******************************/
	 /* solve differential equation */
	 /*******************************/
    printf("Grav. const: %lf\n", G_const);

	 rhs_function_t f;
	 /* here we could do it inplace but this way it is easier to convert to cpn */

	 /* real_t *y = calloc(6*N,sizeof(real_t)); */
	 /* real_t *x = calloc(6*N,sizeof(real_t)); */
	 real_t start_t =start_t_const, end_t=t_end_const, h=h_const;
    real_t t = start_t;


	 f = &nbodyprob;

    print_step(t, output_file,N,particles); 

	 while(t+h < end_t){
		
		  RKstep(f,t,particles,particles,h,N);
		
		  printf("RKstep done: t = %lf\n", t); 

		  print_step(t, output_file,N,particles); 

        t+= h; /* for the non-adapting version */
		  if(t >= end_t){ 
				t = end_t;
		  }
	 }
	 RKstep(f,t,particles,particles,h,N);

	 printf("last RKstep done: t = %lf\n", t);

	 print_step(t, output_file,N,particles); 


	 /************/
	 /* clean up */ 
	 /***********/

	 free(particles);
	 return 0;
}

