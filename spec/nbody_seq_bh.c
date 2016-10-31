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
    for(i=0;i<N;i++){
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
		printf("length = %lf; new length = %lf\n", length, new_length);
            int octant_index = 4*(particle->pos[0] < center.x[0]) + 2*(particle->pos[1] < center.x[1]) + (particle->pos[2] < center.x[2]);
            rvector_t new_center;
            for(j=0;j<3;j++){
printf(" {delta = (%lf - %lf)*%lf}", (real_t)(particle->pos[j] < center.x[j]), 1/2. , new_length);
                new_center.x[j] = center.x[j] + ((real_t)(particle->pos[j] > center.x[j]) - 1/2.) * new_length;
		printf(" %lf ",new_center.x[j]);
            }
printf("\n");

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
		printf("new center for p = [");
		particle_pretty_print(tree->particle);
		printf("] = ");
            for(j=0;j<3;j++){
                new_center.x[j] = center.x[j] + ((real_t)(tree->particle->pos[j] > center.x[j]) - 1/2.) * new_length;
		printf(" %lf ",new_center.x[j]);
            }
printf("\n");

            octree_insert_node(tree->particle,&(tree->children[octant_index]),new_center,new_length);

            /* create a center of mass particle */
            particle_t * center_of_mass = malloc(sizeof(particle_t));
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
	tree->children=NULL;
    }

    /* everything else cleaned now */
    /*free(tree); */ /* this seems to give some trouble now */
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
	 /* res->vel[0] = 0; */
	 /* res->vel[1] = 0; */
	 /* res->vel[2] = 0; */
	 
	 return ret;
}


void read_file(FILE* fp, particle_t *particles){
	 long i = 0;
	 while(read_particle(fp, &(particles[i]) ) == 7) i++;
}


void RKstep ( rhs_function_t f, real_t t, real_t *y, real_t *x, real_t h, long N, particle_t *particles){
  
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

  real_t *ynew = calloc(n*N, sizeof(real_t));
  real_t *sum =  calloc(n*N, sizeof(real_t));
  real_t **K = calloc(m,sizeof(real_t *));
  for(j=0;j<m;j++){
		K[j] = calloc(n*N,sizeof(real_t));
  }/* like: K[m][n*N];  this way K[0] ... K[m-1] are arrays of length N*n */

  /* do{ */ /*this one is for the adaptive version */
	 
    /* init step */ /* unnecesary when using calloc! */
	 /* for(i=0;i<N;i++){ */
		  /* for(l=n*i;l<n*(i+1);l++) */
		  /* 		sum[l] = 0; */

		  /* for(l=n*i;l<n*(i+1);l++) */
        /*     for(j=0;j<m;j++) */
        /*         K[j][l] = 0; */
	 /* } */

	 /* execute the ODE function m times */
	 for(j=0;j<m;j++){

		  for(i=0;i<N;i++)
				for(l=n*i;l<n*(i+1);l++){
					 ynew[l] = y[l];
		  }
		  /* update the arguments to calculate the next K[j] */
		  for(i=0;i<N;i++){

		  
				for(l=n*i;l<n*(i+1);l++){
					 sum[l] = 0;
					 for(k=0;k<j;k++){
						  sum[l]+= Beta[j*m+k]*K[k][l];
					 }
					 /* ynew[l] += h * sum[l] */
					 /* cblas_daxpy(n,h,&sum[l],1,&ynew[l],1); */
					 ynew[l] += h * sum[l];
				}
		  }
		  (*f)(t + Alpha[j] * h, ynew,(real_t *) &K[j][0], N, particles); /* K[j] = (*f)(t + Alpha[j] * h, ynew + (h * sum)); */
		  /*printf("finished. K[%ld][0] = %lf\n",j,K[j][0]); */
	 } 

 
	 for(i=0;i<N;i++){

	   for(l=n*i;l<n*(i+1);l++){
				sum[l] = 0;

		  for(j=0;j<m;j++){
         /* if(i==0) */
         /*     printf("sum[%ld] = %lf*%lf\n",l,Gamma[j],K[j][l]); */
				sum[l]+=Gamma[j]*K[j][l]; //could rewrite in atlas
		  }

		   ynew[l] = y[l] + h*sum[l]; //atlas?
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
		  for(l=n*i;l<n*(i+1);l++)
				x[l] = ynew[l];
	 
	 /* clean up */

	 for(j=0;j<m;j++){
		  free(K[j]);
	 }
	 free(ynew);
	 free(sum);
	 free(K);

}

/* explicitly not inplace to be able to convert it easier to cpn */
void nbodyprob(real_t t, real_t *y, real_t *x, long N, particle_t *particles){
    long i,j,l;
    real_t norm;
	 real_t soft_eps_sq = soft_eps * soft_eps; 

    const int n = 6; /*components (2*dimensions) */
    /* printf("executing nbody function with t = %lf",t); */
  
    for(i=0;i<N;i++){
        for(l=n*i;l<n*i+n/2;l++) /* the first 3 components */
            x[l] = y[n/2+l];  
        for(l=n*i+n/2;l<n*(i+1);l++) /* the second 3 components (where f will be applied) */
            x[l] = 0; 
        /* if (i==0) printf( "x = %lf. ", x[0]); */
        for(j=0;j<N;j++){
            if(i!=j){

                /* rvector_t diff; */
                /* for(l = 0;l<n/2;l++) diff[l] =  y[n*j+l]-y[n*i+l]; */
                /* norm = pow(cblas_dnrm2(3,diff,1),3);  */
					 norm = 0;
					 for(l = 0; l<n/2;l++){
						  norm += pow(y[n*j+l]-y[n*i+l],2);
						  if(norm < soft_eps_sq) /* minimal distance to avoid numeric errors! */
								norm = soft_eps_sq;
					 }
					 norm = pow(norm,3./2.);

                for(l=0;l<n/2;l++){
                    real_t f_i_j = (particles[j].mass * G_const/norm * (y[n*j+l] - y[n*i+l]));
                    /* if(i==0 && j==1) printf("f_0_1 = %lf * %lf/%lf * (%lf - %lf) =%lf. ",particles[j].mass , G_const,norm , y[n*j+l] , y[n*i+l],f_i_j); */
                    x[n*i + n/2 + l] += f_i_j;
                }
            }
        }
        /* if (i== 0) printf( " (af. int.) x[0] = %lf. ", x[0]); */
    }
}


void print_step(real_t t, real_t *y, char *filename_base, long N, particle_t *particles){
  long i;
  char output_string[strlen(filename_base) + 15];
  strcpy(output_string,filename_base);
  char time[15];
  sprintf(time,"_%010.6lf.out",t);
  strcat(output_string,time);
  struct stat st = {0};

  if (stat("output", &st) == -1) {
		mkdir("output", 0700);
  }
  FILE *output_file = fopen(output_string,"w");

  /* fprintf(output_file,"time t = %lf: \n", t); */
	 for(i=0;i<N;i++){
	   /* fprintf(output_file,"x = %lf, y = %lf, z = %lf, vx = %lf, vy = %lf, vz = %lf \n",y[6*i+0],y[6*i+1],y[6*i+2],y[6*i+3],y[6*i+4],y[6*i+5]); */
        fprintf(output_file,"%lf %lf %lf %lf %lf %lf %lf \n",y[6*i+0],y[6*i+1],y[6*i+2],y[6*i+3],y[6*i+4],y[6*i+5],particles[i].mass);
	 }
  fclose(output_file);
}

int main(){


	 /****************************/
	 /* read particles from file */
	 /****************************/

	 char *name = "../data/test_octree.tab";
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

	 printf("read a file with %ld particles\n",N);
	 particle_t *particles = calloc(N,sizeof(particle_t));
	 read_file(fp,particles);

    rvector_t zero = { 0., 0., 0.};
    octree_node_t *tree = octree_generate_tree(N, particles, zero,10);
    octree_pretty_print(tree);
    printf("\n");
    octree_free(tree);
	 long i;
	 /* for(i = 0; i < N; i++){ */
	 /* 	  printf("x[%ld]: %lf\n",i,particles[i].pos[0]); */
	 /* } */
	 
	 /*******************************/
	 /* solve differential equation */
	 /*******************************/
    printf("Grav. const: %lf\n", G_const);

	 rhs_function_t f;
	 /* here we could do it inplace but this way it is easier to convert to cpn */

	 real_t *y = calloc(6*N,sizeof(real_t));
	 /* real_t *x = calloc(6*N,sizeof(real_t)); */
	 real_t start_t =0, end_t=10, h=0.1;
    real_t t = start_t;
	 int l;

    /* for(i=0;i<6*N;i++) x[i] = 0; */
    for(i=0;i<N;i++){ 
        /* printf("particles[%ld].x = %lf",i,particles[i].pos[0]); */
		  for(l=0;l<3;l++){
            y[6*i+l] = particles[i].pos[l];
            y[6*i+3+l] = particles[i].vel[l];
        }
    }


	 f = &nbodyprob;

    print_step(t, y, output_file,N,particles); 

	 while(t+h < end_t){
		
		  RKstep(f,t,y,y,h,N,particles);
		
		  printf("RKstep done: t = %lf\n", t); 

		  print_step(t, y, output_file,N,particles); 

        t+= h; /* for the non-adapting version */
		  if(t >= end_t){ 
				t = end_t;
		  }
	 }
	 RKstep(f,t,y,y,h,N,particles);

	 printf("last RKstep done: t = %lf\n", t);

	 print_step(t, y,output_file,N,particles); 


	 /************/
	 /* clean up */ 
	 /***********/

	 free(particles);
	 return 0;
}

