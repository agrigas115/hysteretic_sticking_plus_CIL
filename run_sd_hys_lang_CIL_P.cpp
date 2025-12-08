/* Run NVT */
/* Alex Grigas */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h> 

#include "voro++.hh"
#include <vector>
#include <iostream>
#include <cmath>

using namespace voro;
using namespace std;

#define pi 3.14159265358979323846

// Cell-cell interaction
// Hysteretic sticky spring
double compute_slj(double coords[][2], int nonbonded_array1[], int nonbonded_array2[], 
					double total_force[][2], int num_neighs, int vlist[], double Lx, double Ly, 
					double vstress[][2], double sigma_ij_array[], int num_atoms, int adhesion_array[], int* z){
	double distance, mag;
	double V = 0;
	int i, index;
	double x1, x2, y1, y2, delta_x, delta_y;
	double sigma_ij;
	*z = 0;
	double delta_x_norm, delta_y_norm, inv_dist;
	double Fx, Fy;
	// Check Verlet neighbors
	for (i = 0; i < num_neighs; i++){
		index = vlist[i];
		x1 = coords[nonbonded_array1[index]][0];
		y1 = coords[nonbonded_array1[index]][1];
		x2 = coords[nonbonded_array2[index]][0];
		y2 = coords[nonbonded_array2[index]][1];

		delta_x = x1-x2;
		delta_x -= Lx * nearbyint(delta_x / Lx);
		
		delta_y = y1-y2;
		delta_y -= Ly * nearbyint(delta_y / Ly);

		distance = sqrt((delta_x*delta_x) + (delta_y*delta_y));

		// Sum of radii
		sigma_ij = sigma_ij_array[index];
		
		// If bonded
		if (adhesion_array[index] == 1){
			// Linear spring magnitude
			mag = (sigma_ij - distance);
			*z += 1;
			if (distance < sigma_ij){
				mag *= 100;
			}
		}
		// If not bonded, check for overlap
		else if (adhesion_array[index] == 0){
			// If overlapping
			if (distance <= sigma_ij){
				mag = (sigma_ij - distance);
				mag *= 100;
				// Form double sided spring
				adhesion_array[index] = 1;
			}
			else{
				mag = 0.0;
			}
		}

		// Get unit vector
		inv_dist = 1.0 / distance;
		delta_x_norm = delta_x * inv_dist;
		delta_y_norm = delta_y * inv_dist;

		// Scale by magnitude
		Fx = delta_x_norm * mag;
		Fy = delta_y_norm * mag;
		
		// Add to the total force
		total_force[nonbonded_array1[index]][0] += Fx;
		total_force[nonbonded_array1[index]][1] += Fy;
		total_force[nonbonded_array2[index]][0] += -Fx;
		total_force[nonbonded_array2[index]][1] += -Fy;

		// Compute Vstress
		vstress[0][0] += Fx*delta_x;
		vstress[1][1] += Fy*delta_y;

		vstress[0][1] += Fx*delta_y;
		vstress[1][0] += Fy*delta_x;

	}
	return V;
}

// CIL-P activity
void update_bonds(double coords[][2], int adhesion_array[], double p_on, double p_off, int num_neighs, 
	int vlist[], double Lx, double Ly, int nonbonded_array1[], int nonbonded_array2[], 
	int num_atoms, int adhesion_idx_array[][25], double adhesion_theta_array[][25], int count_array[], 
	double theta_cutoff, int new_adhesion_array[], int neighbor_array[256][256], double r0_array[], double ratchet_rate, double sigma_ij_array[], double contraction_array[]){

	int i, j, index, count, idx, idx_i, idx_j, count_i, count_j;
	double rand_unif;
	double x1, x2, y1, y2, delta_x, delta_y, distance;
	double theta;
	
	double delta_theta;
	double theta_low, theta_high;

	double sep_theta;
	int angle_1_sat = 0;
	int angle_2_sat = 0;

	int n;
    bool swapped;
    double temp;

    // Check unbinding
    for (i = 0; i < num_neighs; i++){
		index = vlist[i];
		if (adhesion_array[index] == 1){

			rand_unif = drand48();
			if (rand_unif < p_off){
				adhesion_array[index] = 0;
				contraction_array[index] = 0;
			}
		}
	}

	// Set up adhesion pairs
	for (idx = 0; idx < num_neighs; idx++){
		index = vlist[idx];
		if (adhesion_array[index] == 1){
			i = nonbonded_array1[index];
			j = nonbonded_array2[index];
			count = count_array[i];
			adhesion_idx_array[i][count] = j;
			count_array[i] += 1;

			i = nonbonded_array2[index];
			j = nonbonded_array1[index];
			count = count_array[i];
			adhesion_idx_array[i][count] = j;
			count_array[i] += 1;
		}
	}

	// Get pairwise sorted angle lists
	for (i = 0; i < num_atoms; i++){
		count = count_array[i];
		for (idx=0; idx<count; idx++){
			j = adhesion_idx_array[i][idx];

			x1 = coords[i][0];
			y1 = coords[i][1];
			x2 = coords[j][0];
			y2 = coords[j][1];

			delta_x = x2-x1;
			delta_x -= Lx * nearbyint(delta_x / Lx);
			
			delta_y = y2-y1;
			delta_y -= Ly * nearbyint(delta_y / Ly);

			theta = atan2(delta_y, delta_x);

			adhesion_theta_array[i][idx] = theta;
		}

		n = count;
	  	// Bubble sort row
	    for (int idx = 0; idx < n - 1; idx++) {
	        swapped = false;
	        for (int j = 0; j < n - idx - 1; j++) {
	            if (adhesion_theta_array[i][j] > adhesion_theta_array[i][j + 1]) {
	            	// Swap
	                temp = adhesion_theta_array[i][j];
	                adhesion_theta_array[i][j] = adhesion_theta_array[i][j + 1];
	                adhesion_theta_array[i][j + 1] = temp;
	                swapped = true;
	            }
	        }
	      
	        // If no two elements were swapped, then break
	        if (!swapped)
	            break;
	    }
	    
	}

	count = 0;
	for (i=0;i<num_neighs;i++){
		index = vlist[i];
		idx_i = nonbonded_array1[index];
		idx_j = nonbonded_array2[index];
		if (adhesion_array[count] == 0 && neighbor_array[idx_i][idx_j] == 1){

			count_i = count_array[idx_i];
			count_j = count_array[idx_j];

			angle_1_sat = 0;
			angle_2_sat = 0;

			x1 = coords[idx_i][0];
			y1 = coords[idx_i][1];
			x2 = coords[idx_j][0];
			y2 = coords[idx_j][1];

			delta_x = x2-x1;
			delta_x -= Lx * nearbyint(delta_x / Lx);
			
			delta_y = y2-y1;
			delta_y -= Ly * nearbyint(delta_y / Ly);
			distance = sqrt((delta_x*delta_x) + (delta_y*delta_y));

			sep_theta = atan2(delta_y, delta_x);

			// Search if rij is between a large enough angle
			// If only one adhesion, accept
			if (count_i > 1){

				for (idx=0;idx<count_i-1;idx++){

					theta_low = adhesion_theta_array[idx_i][idx];
					theta_high = adhesion_theta_array[idx_i][idx+1];
					delta_theta = theta_high - theta_low;

					if (delta_theta > theta_cutoff){
						if (theta_low < sep_theta && sep_theta < theta_high){
							angle_1_sat = 1;
						}
					}
					
				}
				
				theta_high = adhesion_theta_array[idx_i][count_i-1];
				theta_low = adhesion_theta_array[idx_i][0];

				delta_theta = 2*pi - (theta_high - theta_low);
				if (delta_theta > theta_cutoff){
					if (theta_low > sep_theta || sep_theta > theta_high){
						angle_1_sat = 1;
					}
				}
			
			}
			else{
				angle_1_sat = 1;
			}
			

			// Search other cell

			sep_theta = atan2(-delta_y, -delta_x);

			// Search if rij is between a large enough angle
			
			if (count_j > 1){
				for (idx=0;idx<count_j-1;idx++){
					theta_low = adhesion_theta_array[idx_j][idx];
					theta_high = adhesion_theta_array[idx_j][idx+1];
					delta_theta = theta_high - theta_low;

					if (delta_theta > theta_cutoff){
						if (theta_low < sep_theta && sep_theta < theta_high){
							angle_2_sat = 1;
						}
					}
					
				}
				theta_high = adhesion_theta_array[idx_j][count_j-1];
				theta_low = adhesion_theta_array[idx_j][0];

				delta_theta = 2*pi - (theta_high - theta_low);
				if (delta_theta > theta_cutoff){
					if (theta_low > sep_theta || sep_theta > theta_high){
						angle_2_sat = 1;

					}
				}
			}
			else{
				angle_2_sat = 1;
			}
			
			if ((angle_1_sat == 1 && angle_2_sat == 1) && new_adhesion_array[idx_i] == 0 && new_adhesion_array[idx_j] == 0){
				rand_unif = drand48();
				if (rand_unif < p_on){
					adhesion_array[index] = 1;
					new_adhesion_array[idx_i] = 1;
					new_adhesion_array[idx_j] = 1;
					r0_array[index] = distance;
				}
			}
		}
		count += 1;
	}
	
	
	return;
}


int generate_vlist(double coords[][2], int nonbonded_array1[], int nonbonded_array2[],
					int num_nonbonded, int vlist[], double r_l, double Lx, double Ly, 
					double sigma_ij_array[], int adhesion_array[], int num_atoms){
	int i, idx_i, idx_j;
	int num_neighs = 0;
	double x1, x2, y1, y2, delta_x, delta_y;
	double distance;
	double cutoff;
	for (i=0; i<num_nonbonded; i++){
		if (adhesion_array[i] == 1){
			vlist[num_neighs] = i;
			num_neighs += 1;
		}
		else{
			idx_i = nonbonded_array1[i];
			idx_j = nonbonded_array2[i];
			x1 = coords[idx_i][0];
			y1 = coords[idx_i][1];
			x2 = coords[idx_j][0];
			y2 = coords[idx_j][1];

			cutoff = sigma_ij_array[i] * r_l;

			delta_x = x1-x2;
			delta_x -= Lx * nearbyint(delta_x / Lx);
			
			delta_y = y1-y2;
			delta_y -= Ly * nearbyint(delta_y / Ly);

			distance = sqrt((delta_x*delta_x) + (delta_y*delta_y));

			if (distance < cutoff){
				vlist[num_neighs] = i;
				num_neighs += 1;
			}
		}
		
	}
	return num_neighs;
}



int generate_vlist_2(double coords[][2], int nonbonded_array1[], int nonbonded_array2[],
					int num_nonbonded, int vlist[], double r_l, double Lx, double Ly, 
					double sigma_ij_array[], int adhesion_array[], int num_atoms, int neighbor_array[256][256]){
	int i, idx_i, idx_j;
	int num_neighs = 0;
	double x1, x2, y1, y2, delta_x, delta_y;
	double distance;
	double cutoff;
	for (i=0; i<num_nonbonded; i++){
		idx_i = nonbonded_array1[i];
		idx_j = nonbonded_array2[i];

		if (adhesion_array[i] == 1 || neighbor_array[idx_i][idx_j]){
			vlist[num_neighs] = i;
			num_neighs += 1;
		}
		else{
			x1 = coords[idx_i][0];
			y1 = coords[idx_i][1];
			x2 = coords[idx_j][0];
			y2 = coords[idx_j][1];

			cutoff = sigma_ij_array[i] * r_l;

			delta_x = x1-x2;
			delta_x -= Lx * nearbyint(delta_x / Lx);
			
			delta_y = y1-y2;
			delta_y -= Ly * nearbyint(delta_y / Ly);

			distance = sqrt((delta_x*delta_x) + (delta_y*delta_y));

			if (distance < cutoff){
				vlist[num_neighs] = i;
				num_neighs += 1;
			}
		}
		
	}
	return num_neighs;
}


void generate_gaussian(int half_sample_size, double gaussian_array[][2]){
	int i, dim;
	double rand_unif1, rand_unif2;
	double rand_gauss1, rand_gauss2;
	double twopi = 6.28318530718;
	double A, B;
	for (dim = 0; dim < 2; dim++){
		for (i = 0; i < half_sample_size; i++){
			rand_unif1 = drand48();;
			rand_unif2 = drand48();;

			A = sqrt(-2.*log(rand_unif1));
			B = twopi * rand_unif2;
			rand_gauss1 = A*cos(B);
			rand_gauss2 = A*sin(B);

			gaussian_array[i*2][dim] = rand_gauss1;
			gaussian_array[(i*2)+1][dim] = rand_gauss2;
		}
	}
}

int main(int argc, char *argv[]) {
	// argv[1]: int for sequence
	// argv[2]: run number
	// argv[3]: pon - probability of making a new adhesion
	// argv[4]: poff - probability of cutting an exisiting adhesion
	// argv[5]: phi - packing fraction
	// argv[6]: T - temperature in Langevin thermostat
	// argv[7]: gamma - friction in Langevin thermostat
	// argv[8]: theta - cutoff in adhesion angle to determine if a new adhesion is possible
	// argv[9]: rr - ratchet rate of pulling in adhesion rest length

	srand48(atoi(argv[1])+atoi(argv[2])); // seed RNG based on run number
	char buf[0x100]; // buffering strings

	//printf("Welcome!\n");
	double avg_x, avg_y; // centering
	//clock_t time_1, time_2; // timing
	int i, j, k; // indices
	int num_neighs; // verlet-list
	FILE *fp; // opening files
	double r_l = 1.1; // Multiplier for verlet-list skin
	double dt = 0.1; // Time step
	double dt_lang = dt;

	int z;

	double vstress[2][2] = {0};
	double vstress_t[2][2] = {0};

	//double initial_temp, temp;

	int num_atoms = 256;  
    int half_num_atoms = num_atoms / 2;

	// Initialize Coords //
	double coords[num_atoms][2], half_coords[num_atoms][2], displacement_array[num_atoms][2];
	double total_force[num_atoms][2];
	double velocs[num_atoms][2], half_velocs[num_atoms][2], v_prime[num_atoms][2];
	double sigma_i_array[num_atoms];

	// Load particle sizes //
	snprintf(buf, sizeof(buf), "seqs/seq_%s.txt",argv[1]);
	fp = fopen(buf, "r");
	for (i=0; i<num_atoms; i++){
		fscanf(fp," %lf",&sigma_i_array[i]);
	}
	fclose(fp);

	
	double p_on = atof(argv[3]);
	double p_off = atof(argv[4]);
	
	double target_phi = atof(argv[5]);
	double target_temp = atof(argv[6]);

	double gamma = atof(argv[7]);

	double theta_cutoff = atof(argv[8]);
	double ratchet_rate = atof(argv[9]);
	ratchet_rate /= 1000;

	int adhesion_idx_array[256][25];
	double adhesion_theta_array[256][25];
	int count_array[256] = {0};
	int new_adhesion_array[256] = {0};
	int neighbor_array[256][256];

	// Find initial Box size //
	double initial_phi = 0.1;
	double vol = 0;
	for (i=0; i<num_atoms; i++){
		vol += pi*sigma_i_array[i]*sigma_i_array[i];
	}
	double box_size = sqrt(vol / initial_phi);

	// Assume unit mass


//////////

	// Loading Nonbonded Interactions //
	// Count number of nonboned pairs
	int num_nonbonded = 0; 
	for (i=0; i<num_atoms; i++){
		for(j=i+1; j<num_atoms; j++){
			num_nonbonded += 1;
		}
	}
	// Setting Non-bonded parameters
	double* sigma_ij_array; 
	sigma_ij_array = (double*)malloc(num_nonbonded * sizeof(double));
	double* r0_array; 
	r0_array = (double*)malloc(num_nonbonded * sizeof(double));
	double* tension_array; 
	tension_array = (double*)malloc(num_nonbonded * sizeof(double));
	double* contraction_array; 
	contraction_array = (double*)malloc(num_nonbonded * sizeof(double));

	int count = 0;
	double sigma_ij;
	for (i=0; i<num_atoms; i++){
		for(j=i+1; j<num_atoms; j++){
			sigma_ij = sigma_i_array[i] + sigma_i_array[j];
			sigma_ij_array[count] = sigma_ij;
			r0_array[count] = sigma_ij;
			tension_array[count] = 0;
			contraction_array[count] = 0;
			count += 1;
		}
	}


	// Loading Nonbonded Index Pairs 
	int* adhesion_array; 
	adhesion_array = (int*)malloc(num_nonbonded * sizeof(int));
	int* nonbonded_array1; 
	nonbonded_array1 = (int*)malloc(num_nonbonded * sizeof(int));
	int* nonbonded_array2; 
	nonbonded_array2 = (int*)malloc(num_nonbonded * sizeof(int));
	
	count = 0;
	for (i=0; i<num_atoms; i++){
		for(j=i+1; j<num_atoms; j++){
			nonbonded_array1[count] = i;
			nonbonded_array2[count] = j;
			adhesion_array[count] = 0;
			count += 1;
		}
	}

	// Initialize vlist //
	int* vlist; 
	vlist = (int*)malloc(num_nonbonded * sizeof(int));

	int restart;
	double Lx, Ly;
	double gaussian_array[num_atoms][2];
	int NVT_count = -1;

	// Check if a checkpoint file exists
	snprintf(buf, sizeof(buf), "chk/coords_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.chk",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
	fp = fopen(buf, "r");
	if (fp != NULL){
		restart = 1;
		// Load Step num, box size, and coords
		if (fscanf(fp," %d",&NVT_count)==0){
	    	printf("Problem\n");
	    }
	    if (fscanf(fp," %lf",&Lx)==0){
	    	printf("Problem\n");
	    }
	    if (fscanf(fp," %lf",&Ly)==0){
	    	printf("Problem\n");
	    }
		for (i=0; i<num_atoms; i++){
			for (j=0; j<2; j++){
				if (fscanf(fp," %lf",&coords[i][j])==0){
					printf("Problem\n");
				}
			}
		}
		fclose(fp);

		// Load adhesion array
		snprintf(buf, sizeof(buf), "chk/adhesion_array_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.chk",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
		fp = fopen(buf, "r");
		for (i=0; i<num_nonbonded; i++){
			if (fscanf(fp," %d",&adhesion_array[i])==0){
				printf("Problem\n");
			}
		}

		// Load r0 array
		snprintf(buf, sizeof(buf), "chk/r0_array_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.chk",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
		fp = fopen(buf, "r");
		for (i=0; i<num_nonbonded; i++){
			if (fscanf(fp," %lf",&r0_array[i])==0){
				printf("Problem\n");
			}
		}
	}
	else{
		restart = 0;
		// If not from check point, delete old data files
		snprintf(buf, sizeof(buf), "vstress_data/vstress_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
		remove(buf);
		snprintf(buf, sizeof(buf), "vstress_data/tension_contraction_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
		remove(buf);
		snprintf(buf, sizeof(buf), "bonded_data/bonded_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
		remove(buf);
		snprintf(buf, sizeof(buf), "bonded_data/tension_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
		remove(buf);
		snprintf(buf, sizeof(buf), "traj_data/traj_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
		remove(buf);

		// Initialize random dilute packing
		box_size = sqrt(vol / initial_phi);
		Lx = box_size;
		Ly = box_size;

		/* Getting Initial Coords - RW */
		
		memset(coords, 0, sizeof(coords));
		memset(gaussian_array, 0, sizeof(gaussian_array));
		double delta_x, delta_y;
		//generate_gaussian(half_sample_size, gaussian_array);
		double random_vec_x = 0, random_vec_y = 0;
		int overlap_check = 0, overlap_count = 0;
		double distance, x1, x2, y1, y2;
		for (i=0; i<num_atoms; i++){
			overlap_check = 0;
			// Normalize and scale
			while (overlap_check == 0){
				overlap_count = 0;
				//generate_gaussian(half_sample_size, gaussian_array);
				random_vec_x = drand48();
				random_vec_y = drand48();
				//if (i == 0){
					//coords[i][0] = 0.;
					//coords[i][1] = 0.;
					//overlap_check = 1;
				//}
				//else{

				coords[i][0] = (box_size * random_vec_x) - (box_size/2.);
				coords[i][1] = (box_size * random_vec_y) - (box_size/2.);

				//Check for overlaps
				for (j=0; j<i; j++){
					
					x1 = coords[i][0];
					y1 = coords[i][1];
					x2 = coords[j][0];
					y2 = coords[j][1];

					delta_x = x1-x2;
					delta_x -= box_size * nearbyint(delta_x / box_size);
					
					delta_y = y1-y2;
					delta_y -= box_size * nearbyint(delta_y / box_size);

					sigma_ij = sigma_i_array[i] + sigma_i_array[j];

					distance = sqrt(delta_x*delta_x + delta_y*delta_y);
					if (distance/sigma_ij < 1){
						overlap_check = 0;
						overlap_count += 1;
					}
					
				}
				if (overlap_count == 0){
					overlap_check = 1;
				}
				//}
			}	
		}
	}

	generate_gaussian(half_num_atoms, gaussian_array);
	for (i=0; i<num_atoms; i++){
		for (j=0; j<2; j++){
			total_force[i][j] = 0.;
			displacement_array[i][j] = 0.;
			velocs[i][j] = 0.;
		}
	}


	// Begin Velocity Verlet //
	memset(displacement_array, 0, sizeof(displacement_array));
	memset(total_force, 0, sizeof(total_force));
	memset(vstress, 0, sizeof(vstress));

	/*
	total_kinetic_energy = 0;
	for (j = 0; j < num_atoms; j++){
		v = sqrt( (velocs[j][0]*velocs[j][0]) + (velocs[j][1]*velocs[j][1]) );
		total_kinetic_energy += 0.5 * v * v;
	}
	temp = total_kinetic_energy / num_atoms;

	scaling = sqrt(target_temp/temp);
	for (j = 0; j < num_atoms; j++){
		for (k = 0; k < 2; k++){
			velocs[j][k] *= scaling;
		}
	}
	*/

	// NVE //
	double disp = 0, disp_max1 = 0, disp_max2 = 0;
	double disp_x, disp_y;
	double disp_list[num_atoms];
	double min_rc = 0.5;
	int split = 1000;

	num_neighs = generate_vlist(coords, nonbonded_array1, nonbonded_array2, num_nonbonded, vlist, r_l, Lx, Ly, sigma_ij_array, adhesion_array, num_atoms);
	compute_slj(coords, nonbonded_array1, nonbonded_array2, total_force, num_neighs, vlist, Lx, Ly, vstress, sigma_ij_array, num_atoms, adhesion_array, &z,r0_array, ratchet_rate, tension_array, contraction_array);
	i = 0;
	

	// UNCHANGING FIRE VARIABLES
	int npPos = 0;
    int npNeg = 0;
    int npPMin = 0;
    
	double vnorm = 0;
	double fnorm = 0;
	double P = 0;

	double alpha0 = 0.25;
	double finc = 1.01;
	double fdec = 0.25;
	double falpha = 0.99;

	double NNEGMAX = 2000;


    // FIRE Parameters #
    double FIRE_alpha = alpha0;
    double FIRE_dt =  dt;
    
    double dtmax = 0.1;
    double dtmin = 1e-2*FIRE_dt;

	double Ftol = 1e-14;
	double fcheck = 10*Ftol;

	/* Beginning Jamming */
	double pcheck = 0;
	int fireit = 0;
	int it = 0;
	int fireitmin = 100;

	double delta_phi = 1e-3;
	double Lscale, phi, delta_L;

	fcheck = 0;
	for (j = 0; j < num_atoms; j++){
		fcheck += sqrt(total_force[j][0]*total_force[j][0] + total_force[j][1]*total_force[j][1]);
	}
	fcheck = fcheck/num_atoms;
	phi = vol / (Lx*Ly);
	int cycle = 0;
	if (restart == 0){
		while (cycle < 2){
			//printf("phi = %lf\n",phi);
			if (fabs(pcheck) > 1e-12 && cycle == 0){
	    		delta_phi = 1e-3;
	    	}

	    	if (phi >= 0.82){
	    		delta_phi = -1e-3;
	    		cycle += 1;
	    	}
	    	if (phi < target_phi && cycle != 0){
	    		delta_phi = 1e-3;
	    		cycle += 1;
	    	}
			
		    it += 1;

		    Lscale = sqrt(phi/(phi+delta_phi));
		    //printf("Lscale = %lf \n", Lscale);
		    delta_L = box_size - (box_size * Lscale);
		    //printf("delta_L = %lf \n", delta_L);
		    box_size *= Lscale;
		    Lx *= Lscale;
		    Ly *= Lscale;
		    for (i=0; i<num_atoms; i++){
		    	for (j=0; j<2; j++){
		    		coords[i][j] *= 1-(delta_L/box_size);
		    	}
		    }

		    phi = vol / (Lx*Ly);

			num_neighs = generate_vlist(coords, nonbonded_array1, nonbonded_array2, num_nonbonded, vlist, r_l, Lx, Ly, sigma_ij_array, adhesion_array, num_atoms);
			memset(displacement_array, 0, sizeof(displacement_array));

			
			fireit = 0;
			FIRE_dt = dt;
			fcheck = 10*Ftol;
		    while (fcheck > Ftol){
		        fireit += 1;
		        //printf("fireit = %d \n", fireit);
		        // Velocity-verlet step
		        for (i = 0; i<num_atoms; i++){
		        	velocs[i][0] += 0.5 * FIRE_dt * total_force[i][0];
		        	velocs[i][1] += 0.5 * FIRE_dt * total_force[i][1];
		        }   
		        
		        // 1. Calculate P, fnorm and vnorm
		        vnorm = 0;
		        for (i = 0; i < num_atoms; i++){
		        	vnorm += sqrt(velocs[i][0]*velocs[i][0] + velocs[i][1]*velocs[i][1]);
		        }
		        fnorm = 0;
		        for (i = 0; i < num_atoms; i++){
		        	fnorm += sqrt(total_force[i][0]*total_force[i][0] + total_force[i][1]*total_force[i][1]);
		        }
		        P = 0;
		        for (i=0; i < num_atoms; i++){
		        	P += velocs[i][0]*total_force[i][0] + velocs[i][1]*total_force[i][1];
		        }
		        
		        
		        // 2. Adjust simulation based on net motion
		        if (P > 0){
		        	npPos += 1;
		            
		            npNeg = 0;
		            
		            //alphat = alpha;
		           
		            // alter simulation if enough positive steps have been taken
		            if (npPos > fireitmin){
		            	// change time step
		                if (FIRE_dt*finc < dtmax){
		                    FIRE_dt = FIRE_dt*finc;
		                }
		            }
		            // decrease alpha
		            FIRE_alpha = FIRE_alpha*falpha;
		        }
		        else {
		        	// reset positive counter
		            npPos = 0;

		            // increase negative counter
		            npNeg = npNeg + 1;

		            // check for stuck simulation
		            if (npNeg > NNEGMAX){
		                printf("Simulation negative for too long, ending program here \n");
		            }

		            // decrease time step if past initial delay
		            if (fireit > fireitmin){
		                // decrease time step
		                if (FIRE_dt*fdec > dtmin){
		                    FIRE_dt = FIRE_dt*fdec;
		                }

		                // change alpha
		                FIRE_alpha = alpha0;
		            }

		            // take a half step backwards
		            for (i = 0; i<num_atoms; i++){
		            	coords[i][0] -= 0.5*FIRE_dt*velocs[i][0];
		            	coords[i][1] -= 0.5*FIRE_dt*velocs[i][1];
		            }
		            

		            // reset velocities to 0
		            memset(velocs, 0, sizeof(velocs));
		        }
		            
		        
		        // update velocities if forces are acting
		        if (fnorm > 0){
		        	for (i = 0; i<num_atoms; i++){
		        		velocs[i][0] = (1 - FIRE_alpha) * velocs[i][0] + FIRE_alpha*(total_force[i][0]/fnorm)*vnorm;
		        		velocs[i][1] = (1 - FIRE_alpha) * velocs[i][1] + FIRE_alpha*(total_force[i][1]/fnorm)*vnorm;
		        	}
		        }

		        // do first verlet update for vertices (assume unit mass)
		        for (i = 0; i<num_atoms; i++){
		        	coords[i][0] += 0.5*FIRE_dt*velocs[i][0];
		        	coords[i][1] += 0.5*FIRE_dt*velocs[i][1];
		        }


		        // Compute new Force
		        memset(total_force, 0, sizeof(total_force));
		        memset(vstress, 0, sizeof(vstress));

		        // Checking Disp for update on Verlet List 
				for (j=0; j<num_atoms; j++){
					disp_x = displacement_array[j][0]-coords[j][0];
					disp_y = displacement_array[j][1]-coords[j][1];
					disp = sqrt(disp_x*disp_x + disp_y*disp_y);
					disp_list[j] = disp;
				}

				// Find first max 
				disp_max1 = 0;
				for (j=0; j<num_atoms; j++){
					if (disp_list[j] > disp_max1){
						disp_max1 = disp_list[j];
					}
				}

				// Find second max 
				disp_max2 = 0;
				for (j=0; j<num_atoms; j++){
					if (disp_list[j] > disp_max2 && disp_list[j] < disp_max1){
						disp_max2 = disp_list[j];
					}
				}

		        if (disp_max1 + disp_max2 > (min_rc*r_l) - min_rc){
					memcpy(displacement_array, coords, sizeof(displacement_array));
					num_neighs = generate_vlist(coords, nonbonded_array1, nonbonded_array2, num_nonbonded, vlist, r_l, Lx, Ly, sigma_ij_array, adhesion_array, num_atoms);
				}

				compute_slj(coords, nonbonded_array1, nonbonded_array2, total_force, num_neighs, vlist, Lx, Ly, vstress, sigma_ij_array, num_atoms, adhesion_array, &z, r0_array, ratchet_rate, tension_array, contraction_array);

		        // Do second verlet update for velocities
		        for (i = 0; i<num_atoms; i++){
		        	velocs[i][0] += 0.5 * FIRE_dt * total_force[i][0];
		        	velocs[i][1] += 0.5 * FIRE_dt * total_force[i][1];
		        }   
		        // update Kcheck and Fcheck
		        // Calc avg norm of total_force
				fcheck = 0;
				for (j = 0; j < num_atoms; j++){
					fcheck += sqrt(total_force[j][0]*total_force[j][0] + total_force[j][1]*total_force[j][1]);
				}
		    	fcheck = fcheck/num_atoms;

		        if (fcheck < Ftol){
		            npPMin = npPMin + 1;
		        }
		        else{
		            npPMin = 0;
		        }

			}

			pcheck = (vstress[0][0] + vstress[1][1]) / (2*Lx*Ly);
			phi = vol / (Lx*Ly);

		
		}
	}

	double exp_gamma_dt =  exp(-gamma * dt_lang);
	double fluc_coeff = sqrt(1. - exp(-2. * gamma * dt_lang)) * sqrt(target_temp);
	double half_dt_lang = 0.5 * dt_lang;

    phi = vol / (Lx*Ly);
    //printf("phi = %lf\n",phi);

    double x_max = Lx / 2;
    double x_min = -Lx / 2; 
    double y_max = Ly / 2;
    double y_min = -Ly / 2; 

    double volume = Lx * Ly;
	double vx, vy;

	//r_l = 3*Lx/sqrt(num_atoms);
	memcpy(displacement_array, coords, sizeof(displacement_array));
	num_neighs = generate_vlist(coords, nonbonded_array1, nonbonded_array2, num_nonbonded, vlist, r_l, Lx, Ly, sigma_ij_array, adhesion_array, num_atoms);
	while (NVT_count < 1e9){

		
		NVT_count += 1;

		// Update Half Velocities //
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 2; k++){
				half_velocs[j][k] = velocs[j][k] + half_dt_lang * total_force[j][k];
			}
		}

		// Update Half Coordinates //
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 2; k++){
				half_coords[j][k] = coords[j][k] + half_dt_lang * half_velocs[j][k];
			}
		}

		// Calc v_prime //
		generate_gaussian(num_atoms/2, gaussian_array);
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 2; k++){
				v_prime[j][k] = exp_gamma_dt * half_velocs[j][k] + fluc_coeff * gaussian_array[j][k];
			}
		}

		// Update Coordinates //
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 2; k++){
				coords[j][k] = half_coords[j][k] + half_dt_lang * v_prime[j][k];
			}
		}
		// Compute new Force
        memset(total_force, 0, sizeof(total_force));
        memset(vstress, 0, sizeof(vstress));

        // Checking Disp for update on Verlet List 
		for (j=0; j<num_atoms; j++){
			disp_x = displacement_array[j][0]-coords[j][0];
			disp_y = displacement_array[j][1]-coords[j][1];
			disp = sqrt(disp_x*disp_x + disp_y*disp_y);
			disp_list[j] = disp;
		}

		// Find first max 
		disp_max1 = 0;
		for (j=0; j<num_atoms; j++){
			if (disp_list[j] > disp_max1){
				disp_max1 = disp_list[j];
			}
		}

		// Find second max 
		disp_max2 = 0;
		for (j=0; j<num_atoms; j++){
			if (disp_list[j] > disp_max2 && disp_list[j] < disp_max1){
				disp_max2 = disp_list[j];
			}
		}

        if (disp_max1 + disp_max2 > (min_rc*r_l) - min_rc){
			memcpy(displacement_array, coords, sizeof(displacement_array));
			num_neighs = generate_vlist_2(coords, nonbonded_array1, nonbonded_array2, num_nonbonded, vlist, r_l, Lx, Ly, sigma_ij_array, adhesion_array, num_atoms, neighbor_array);
		}

		compute_slj(coords, nonbonded_array1, nonbonded_array2, total_force, num_neighs, vlist, Lx, Ly, vstress, sigma_ij_array, num_atoms, adhesion_array, &z, r0_array, ratchet_rate, tension_array, contraction_array);

		// Update Velocities //
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 2; k++){
				velocs[j][k] = v_prime[j][k] + half_dt_lang * total_force[j][k];
			}
		} 

		if (NVT_count % split == 0){

			memset(neighbor_array, 0, sizeof(neighbor_array));

		    // Run Voro++
			int n_x = 5, n_y = 5, n_z = 5; // Grid subdivisions
		    bool periodic_x = true, periodic_y = true, periodic_z = false; // Periodicity

		    // Create a container
		    container_poly con(x_min, x_max, y_min, y_max, -0.5, 0.5, 
                  n_x, n_y, n_z, periodic_x, periodic_y, periodic_z, 8);

			for (i=0;i<num_atoms;i++){
				con.put(i,coords[i][0],coords[i][1],0,sigma_i_array[i]);
			}

			voronoicell_neighbor c;
		    c_loop_all cl(con);
		    if (cl.start()) do {
		        int id = cl.pid();

		        if (con.compute_cell(c, cl)) { // Compute Voronoi cell
		            std::vector<int> neighbors;
		            c.neighbors(neighbors); // Get neighbor IDs
		            
		            // Print neighbors
		            
		            for ( i = 0; i < neighbors.size(); i++) {
    					int nid = neighbors[i];
    					if (nid >= 0){
    						neighbor_array[id][nid] = 1;
    					}
		                
		            }
		            
		        }
		    } while (cl.inc());

			memset(adhesion_idx_array, 0, sizeof(adhesion_idx_array));
			memset(adhesion_theta_array, 0, sizeof(adhesion_theta_array));
			memset(count_array, 0, sizeof(count_array));
			memset(new_adhesion_array, 0, sizeof(new_adhesion_array));
			update_bonds(coords, adhesion_array, p_on, p_off, num_neighs, vlist, Lx, Ly, nonbonded_array1, nonbonded_array2, num_atoms, adhesion_idx_array, adhesion_theta_array, count_array, theta_cutoff, new_adhesion_array, neighbor_array, r0_array, ratchet_rate, sigma_ij_array, contraction_array);
		
			memcpy(displacement_array, coords, sizeof(displacement_array));
			num_neighs = generate_vlist_2(coords, nonbonded_array1, nonbonded_array2, num_nonbonded, vlist, r_l, Lx, Ly, sigma_ij_array, adhesion_array, num_atoms, neighbor_array);

		}


		if (NVT_count % 100000 == 0){
			//printf("NVT_count = %d\n",NVT_count);

			// Get Temperature stress tensor //
			memset(vstress_t, 0, sizeof(vstress));

			// Average Velocs //
			avg_x = 0;
			avg_y = 0;
			for (j=0; j<num_atoms; j++){
				avg_x += velocs[j][0];
				avg_y += velocs[j][1];
			}
			avg_x /= num_atoms;
			avg_y /= num_atoms;

			// Add Veloc comp. to vstress //
			for (j = 0; j < num_atoms; j++){

				vx = velocs[j][0]-avg_x;
				vy = velocs[j][1]-avg_y;

				vstress_t[0][0] += vx * vx;
				vstress_t[1][1] += vy * vy;

				vstress_t[0][1] += vx * vy;
				vstress_t[1][0] += vy * vx;

			}
			
			// Normalized stress tensors by system volume //
			vstress[0][0] /= volume;
			vstress[1][1] /= volume;
			vstress[0][1] /= volume;
			vstress[1][0] /= volume;

			vstress_t[0][0] /= volume;
			vstress_t[1][1] /= volume;
			vstress_t[0][1] /= volume;
			vstress_t[1][0] /= volume;

			snprintf(buf, sizeof(buf), "vstress_data/vstress_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
			FILE *out_file_vstress = fopen(buf, "a+");
			fprintf(out_file_vstress, "%.16g\t%.16g\t%.16g\t%.16g\t%.16g\t%.16g\t%d\n",vstress[0][0],vstress[1][1],vstress[0][1],vstress_t[0][0],vstress_t[1][1],vstress_t[0][1],z);
			fclose(out_file_vstress);

			/*
			snprintf(buf, sizeof(buf), "vstress_data/tension_contraction_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
			FILE *out_file_contract = fopen(buf, "a+");
			for (j = 0; j < num_nonbonded; j++){
				if (adhesion_array[j] == 1){
					fprintf(out_file_contract, "%.16g\t%.16g\n",tension_array[j], contraction_array[j]);
				}
			}
			fclose(out_file_contract);
			*/

			snprintf(buf, sizeof(buf), "bonded_data/bonded_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
			FILE *out_file_bonded = fopen(buf, "a+");
			for (j = 0; j < num_nonbonded; j++){
				if (adhesion_array[j] == 1){
					fprintf(out_file_bonded, "%d\t%d\t",nonbonded_array1[j], nonbonded_array2[j]);
				}
			}
			fprintf(out_file_bonded, "\n");
			fclose(out_file_bonded);

			/*
			snprintf(buf, sizeof(buf), "bonded_data/tension_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
			FILE *out_file_tension = fopen(buf, "a+");
			for (j = 0; j < num_nonbonded; j++){
				if (adhesion_array[j] == 1){
					fprintf(out_file_tension, "%.5g\t%.5g\t",tension_array[j], r0_array[j]);
				}
			}
			fprintf(out_file_tension, "\n");
			fclose(out_file_tension);
			*/
			
			snprintf(buf, sizeof(buf), "traj_data/traj_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
			FILE* out_file_coords_final = fopen(buf, "a+"); // write only 
			fprintf(out_file_coords_final, "%.10g\t%.10g\n", Lx, Ly);
			for (j = 0; j < num_atoms; j++){
				fprintf(out_file_coords_final, "%.10g\t%.10g\n", coords[j][0], coords[j][1]); // write to file
			}
			fclose(out_file_coords_final);
			
			
			snprintf(buf, sizeof(buf), "chk/coords_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.chk",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
			FILE* out_file_coords_chk = fopen(buf, "w"); // write only 
			fprintf(out_file_coords_chk, "%d\t%.16g\t%.16g\n", NVT_count, Lx, Ly);
			for (j = 0; j < num_atoms; j++){
				fprintf(out_file_coords_chk, "%.16g\t%.16g\n", coords[j][0], coords[j][1]); // write to file
			}
			fclose(out_file_coords_chk);

			snprintf(buf, sizeof(buf), "chk/adhesion_array_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.chk",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
			FILE* out_file_adhesion_chk = fopen(buf, "w"); // write only
			for (j = 0; j < num_nonbonded; j++){
				fprintf(out_file_adhesion_chk, "%d\n", adhesion_array[j]);
			}
			fclose(out_file_adhesion_chk);

			snprintf(buf, sizeof(buf), "chk/r0_array_seq_%s_%s_pon_%s_poff_%s_phi%s_T_%s_gamma_%s_theta_cutoff_%s_rr_%s.chk",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
			FILE* out_file_r0_chk = fopen(buf, "w"); // write only
			for (j = 0; j < num_nonbonded; j++){
				fprintf(out_file_r0_chk, "%.16g\n", r0_array[j]);
			}
			fclose(out_file_r0_chk);
			
		}
	}

}

