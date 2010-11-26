/**
 *******************************************************************************
 * Copyright (C) 2010 Nicolas Martinelli, nicolas.martinelli@gmail.com         *
 * Adapted from routines done by Yoann Olivier.                                *
 *                                                                             *
 * This program is free software: you can redistribute it and/or modify        *
 * it under the terms of the GNU General Public License as published by        *
 * the Free Software Foundation, either version 3 of the License, or           *
 * (at your option) any later version.                                         *
 *                                                                             *
 * This program is distributed in the hope that it will be useful,             *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of              *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
 * GNU General Public License for more details.                                *
 *                                                                             *
 * You should have received a copy of the GNU General Public License           *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>        *
 *******************************************************************************
 */

// To compile on lucky : 
// g++ -Wall MC_BKL_layer.cpp /home/nmartine/lib/boost/lib/libboost_thread.so -O3 -lm -I /home/nmartine/lib/boost/include/ -Wl,-rpath,/home/nmartine/lib/boost/lib -o MC_BKL_layer
//#include "/home/nmartine/lib/boost_1_43_0/boost/thread/thread.hpp"
//#include <boost/thread.hpp>

#define _INSIDE_MC_BKL

// C++ libraries
#include <iostream> //Entrées-sorties standard
#include <fstream> //Entrées-sorties sur fichiers
#include <string> //Chaines de caracteres
#include <limits> //Pour aller à la fin d'une ligne, par exemple
#include <iomanip> //Manipulation des sorties, par exemple format scientifique
#include <sstream> //Pour les conversion entre types de variables
#include <vector>
#include <algorithm>

// C libraries
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <omp.h>

// Local libraries
#include "constants.h"
#include "variables.h"
#include "clear.h"
#include "coordinates.h"
#include "mathplus.h"
#include "read.h"
#include "buildgrid.h"
#include "electricfield.h"
#include "distance.h"
#include "deltae.h"
#include "deltav.h"
#include "transferrate.h"
#include "dispatch.h"
#include "printsummary.h"
#include "mcbkl.h"

using namespace std;

// =============================================================================
// ----------------------- Monte-Carlo related functions -----------------------
// =============================================================================

// FRM algorithm
void MC_FRM(string output_folder){
	
	cout << "[INFO] Using the FRM algorithm." << endl;
	
	// Variables for each charge
	vector<int> curr_mol, curr_grid; // Number of the molecule in the mini-grid, and number of this mini-grid
	vector<double> dist, jump; // Distance and number of jumps for each charge
	
	// Variables for a chosen event
	vector<double> event_k; // Transfer rate
	vector<int> event_charge, event_mol_index, event_neigh_num, event_neigh_index; // Initial molecule and neighbor
	
	// Total time and distance
	double total_time_try, total_dist_try;
	
	// Average mobility for each frame
	vector<double> mu_frame;
	
	double uF_x, uF_y, uF_z;
	uF_x = F_x/F_norm;
	uF_y = F_y/F_norm;
	uF_z = F_z/F_norm;
	
	// Check if it's possible to write output files
	stringstream OUT_SIMU, OUT_ERROR, T, P;
  
	// T << theta_deg;
	// P << phi_deg;
	// OUT_SIMU << output_folder << "/simu_" << T.str().c_str() << "_" << P.str().c_str() << ".out";
	// OUT_ERROR << output_folder << "/error_" << T.str().c_str() << "_" << P.str().c_str() << ".out";
	
	OUT_SIMU << output_folder << "/simu_" << charge.c_str() << "_" << F_dir.c_str() << ".out";
	OUT_ERROR << output_folder << "/error_" << charge.c_str() << "_" << F_dir.c_str() << ".out";
  
	FILE * pFile;

	pFile = fopen(OUT_SIMU.str().c_str(), "w");
	if (pFile==NULL) {
		int wait = 0; 
		while (wait<10 && pFile==NULL){
			cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) << " seconds to write a file" << endl;
			usleep(10*(wait+1)*(wait+1));
			pFile = fopen(OUT_SIMU.str().c_str(), "w");
			wait++;
		}
		if (wait==10 && pFile==NULL){
			cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
			exit (1);
		}
	}
	fclose(pFile);
	
	// Start the FRM algorithm
	for (int i=0; i<n_frame; i++){
		
		// ---------------------------------------------------------------------
		
		// Print summary for the current try
		pFile=fopen(OUT_SIMU.str().c_str(), "a");
		if (pFile==NULL) {
			int wait = 0; 
			while (wait<10 && pFile==NULL){
				cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) << " seconds to write a file" << endl;
				usleep(10*(wait+1)*(wait+1));
				pFile=fopen(OUT_SIMU.str().c_str(), "a");
				wait++;
			}
			if (wait==10 && pFile==NULL){
				cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
				exit (1);
			}
		}
		fprintf(pFile,"===============================================================================\n");
		fprintf(pFile,"-------------------------------------------------------------------------------\n");
		fclose(pFile);	
		
		// ---------------------------------------------------------------------
		
		// Generate the grid
		grid_occ.clear();
				
		for (int ii=0; ii<n_mol; ii++){
			grid_occ.push_back( vector<bool> ());
			
			for (unsigned int x=0; x<box_neigh_label.size(); x++){
				grid_occ[ii].push_back( false );
			}
		}

		// ---------------------------------------------------------------------

		// Put the charges in the grid
		curr_mol.clear();
		curr_grid.clear();
		int *pos;
		pos = new int[2];
		
		for (unsigned int charge_i=0; charge_i<n_charges; charge_i++){
			if (n_layer == 0)
				Dispatch_Mol_RND(i, grid_occ, pos);
			else
				Dispatch_Mol_RND_layer(i, grid_occ, pos);
			
			curr_mol.push_back(pos[0]);
			curr_grid.push_back(pos[1]);
			grid_occ[pos[0]][pos[1]] = true;
		}

		delete [] pos;
		
		// ---------------------------------------------------------------------

		// Set distances, number of jumps and the travel time to zero
		total_time_try = 0.0;
		total_dist_try = 0.0;

		dist.clear();
		jump.clear(); 
		for (unsigned int charge_i=0; charge_i<n_charges; charge_i++){
			dist.push_back(0.0);
			jump.push_back(0.0);
		}

		bool previous_jump_ok = true;
		bool exit_loop = false;
		
		// ---------------------------------------------------------------------
				
		// Calculates the waiting times for all the charges
		// at the beginning of the simulation
		
		event_k.clear();
		event_charge.clear();
		event_mol_index.clear();
		event_neigh_num.clear(); 
		event_neigh_index.clear();
		
		for (unsigned int charge_i=0; charge_i<curr_mol.size(); charge_i++){
			
			double event_k_tmp = numeric_limits<double>::max();
			int event_charge_tmp = 0;
			int event_mol_index_tmp = 0;
			int event_neigh_num_tmp = 0;
			int event_neigh_index_tmp = 0;
			
			int tmp_mol_index = curr_mol[charge_i];
									
			for (unsigned int jj=0; jj<neigh_label[i][tmp_mol_index].size(); jj++){
				
				// Find index of neighbor
				int tmp_neigh_num = jj;
				int tmp_neigh_index = tmp_mol_index;
				for (int ii=0; ii<n_mol; ii++){
					if (neigh_label[i][tmp_mol_index][jj] == mol_label[ii]){
						tmp_neigh_index = ii;
						break;
					}
				}
				
				// Calculate waiting_time
				double k_tmp = numeric_limits<double>::min();
				if (tmp_neigh_index != tmp_mol_index){
					k_tmp = Marcus_Levich_Jortner_rate_electro(i, tmp_mol_index, tmp_neigh_index, tmp_neigh_num, d_x[i][tmp_mol_index][jj], d_y[i][tmp_mol_index][jj], d_z[i][tmp_mol_index][jj], dE[i][tmp_mol_index][jj], J_H[i][tmp_mol_index][jj], J_L[i][tmp_mol_index][jj], curr_mol, curr_grid, charge_i);
				
					double random_number = Rand_0_1();
					k_tmp = -log(random_number)/(k_tmp);
				
					// Keep the value only if k is larger
					if (k_tmp < event_k_tmp){
						event_k_tmp = k_tmp;
						event_charge_tmp = charge_i;
						event_mol_index_tmp = tmp_mol_index;
						event_neigh_num_tmp = jj;
						event_neigh_index_tmp = tmp_neigh_index;
					}
				}
			}
			
			// Keep only the highest k
			event_k.push_back(event_k_tmp);
			event_charge.push_back(event_charge_tmp);
			event_mol_index.push_back(event_mol_index_tmp);
			event_neigh_num.push_back(event_neigh_num_tmp);
			event_neigh_index.push_back(event_neigh_index_tmp);
		}
		
		// TO DO : SORT THE TABLE ?
		
		bool first_calc = true;
				
		// ---------------------------------------------------------------------

		unsigned int event = 0;
		int charge_previously_selected = 0;
		
		unsigned int charge_try = 0;
		while (charge_try<n_try*n_charges){	
			
			// Calculates the transfer rate for the previously selected molecule
			if (first_calc==false){

				int charge_i = charge_previously_selected;
										
				double event_k_tmp = numeric_limits<double>::max();
				int event_charge_tmp = 0;
				int event_mol_index_tmp = 0;
				int event_neigh_num_tmp = 0;
				int event_neigh_index_tmp = 0;
				
				int tmp_mol_index = curr_mol[charge_i];
										
				for (unsigned int jj=0; jj<neigh_label[i][tmp_mol_index].size(); jj++){
					
					// Find index of neighbor
					int tmp_neigh_num = jj;
					int tmp_neigh_index = tmp_mol_index;
					for (int ii=0; ii<n_mol; ii++){
						if (neigh_label[i][tmp_mol_index][jj] == mol_label[ii]){
							tmp_neigh_index = ii;
							break;
						}
					}
					
					// Calculate waiting_time
					double k_tmp = numeric_limits<double>::max();
					if (tmp_neigh_index != tmp_mol_index){
						k_tmp = Marcus_Levich_Jortner_rate_electro(i, tmp_mol_index, tmp_neigh_index, tmp_neigh_num, d_x[i][tmp_mol_index][jj], d_y[i][tmp_mol_index][jj], d_z[i][tmp_mol_index][jj], dE[i][tmp_mol_index][jj], J_H[i][tmp_mol_index][jj], J_L[i][tmp_mol_index][jj], curr_mol, curr_grid, charge_i);
					
						double random_number = Rand_0_1();
						k_tmp = -log(random_number)/(k_tmp);
				
						// Keep the value only if k is larger
						if (k_tmp < event_k_tmp){
							event_k_tmp = k_tmp;
							event_charge_tmp = charge_i;
							event_mol_index_tmp = tmp_mol_index;
							event_neigh_num_tmp = jj;
							event_neigh_index_tmp = tmp_neigh_index;
						}
					}
				}
				// Keep the highest k
				event_k[charge_previously_selected] = event_k_tmp;
				event_charge[charge_previously_selected] = event_charge_tmp;
				event_mol_index[charge_previously_selected] = event_mol_index_tmp;
				event_neigh_num[charge_previously_selected] = event_neigh_num_tmp;
				event_neigh_index[charge_previously_selected] = event_neigh_index_tmp;
			}
			
			else{
				exit_loop = false;
			}
			
			// Choose the fastest event
			double k_tmp = numeric_limits<double>::max();
			for (unsigned int t=0; t<event_k.size(); t++){
				//cout << event_k[t] << endl;
				if(event_k[t] < k_tmp) {
					k_tmp = event_k[t];
					charge_previously_selected = event_charge[t];
					event = t;
				}
			}

			// =======================================
			//
			// Summary of variables:
			// ---------------------
			//
			// event_charge[event]: index of the charge of the event (= charge_i)
			// event_mol_index[event]: index of the molecule where the chosen charge is (= ii)
			// event_neigh_num[event]: number of the neighbor of the molecule where the chosen charge is (= jj)
			// event_neigh_index[event]: index of the molecule where the chosen charge is going to jump (= ii for the neighbor)
			// curr_mol[event_charge[event]]: index of the molecule where the chosen charge is (= ii)
			// curr_grid[event_charge[event]]: index of the mini-grid where the chosen charge is (= x)
			//
			// =======================================
	
			// Calculate the new position in the grid
			int tmp_curr_mol = 0; // Expected final molecule index
			int tmp_curr_grid = 0; // Expected final mini-grid index
			int tmp_curr_grid_a = 0, tmp_curr_grid_b = 0, tmp_curr_grid_c = 0;
			bool out_of_system = true;
			
			// Find next mini-grid
			if (neigh_jump_vec_a[i][event_mol_index[event]][event_neigh_num[event]] == 0 && neigh_jump_vec_b[i][event_mol_index[event]][event_neigh_num[event]] == 0 && neigh_jump_vec_c[i][event_mol_index[event]][event_neigh_num[event]] == 0){
				tmp_curr_grid = curr_grid[event_charge[event]];
				tmp_curr_grid_a = box_a[tmp_curr_grid];
				tmp_curr_grid_b = box_b[tmp_curr_grid];
				tmp_curr_grid_c = box_c[tmp_curr_grid];
				out_of_system = false;
				
			}
			
			else {
				for (unsigned int xx=0; xx<box_neigh_label[curr_grid[event_charge[event]]].size(); xx++){
					if (neigh_jump_vec_a[i][event_mol_index[event]][event_neigh_num[event]] == box_neigh_a[curr_grid[event_charge[event]]][xx] && neigh_jump_vec_b[i][event_mol_index[event]][event_neigh_num[event]] == box_neigh_b[curr_grid[event_charge[event]]][xx] && neigh_jump_vec_c[i][event_mol_index[event]][event_neigh_num[event]] == box_neigh_c[curr_grid[event_charge[event]]][xx]){
						tmp_curr_grid = box_neigh_label[curr_grid[event_charge[event]]][xx];
						tmp_curr_grid_a = box_a[tmp_curr_grid];
						tmp_curr_grid_b = box_b[tmp_curr_grid];
						tmp_curr_grid_c = box_c[tmp_curr_grid];
						out_of_system = false;
						break;
					}
					
				}
			}
				
			// Find next molecule
			tmp_curr_mol = event_neigh_index[event];
			
			// Do not jump if the charge goes out of the system
			if (out_of_system){
				previous_jump_ok = false;
				exit_loop = true;
				first_calc = false;
				//printf("out of system\n");
			}
			
			// Do not jump if the next molecule is occupied
			else if (grid_occ[tmp_curr_mol][tmp_curr_grid] == true){
				previous_jump_ok = false;
				exit_loop = true;
				first_calc = false;
				//printf("occupied\n");
			}

			// Do not jump if the charge change of layer
			else if (n_layer > 0 && list_layer[i][curr_mol[event_charge[event]]] != list_layer[i][tmp_curr_mol]){
				previous_jump_ok = false;
				exit_loop = true;
				first_calc = false;
			}
			
			else{
				
				// Calculate the total time
				total_time_try += event_k[event];
	
				// Substract the waiting time
				double k_event = event_k[event];
				for (unsigned int t=0; t<event_k.size(); t++){
					event_k[t] = event_k[t] - k_event;
				}

				// Calculate the distance traveled by the charge and the total distance
				double event_dist = (d_x[i][event_mol_index[event]][event_neigh_num[event]] * uF_x + d_y[i][event_mol_index[event]][event_neigh_num[event]] * uF_y + d_z[i][event_mol_index[event]][event_neigh_num[event]] * uF_z)*1E-8;
				dist[event_charge[event]] += event_dist;
				//total_dist_try += event_dist;
				total_dist_try += event_dist/double(curr_mol.size());
				// Note:
				// We divide here by the number of charges in the system, 
				// so there is no problem afterwards if the number of charges
				// in the system change. Indeed:
				// <x> = (1/n) * sum_i x_i = sum_i (x_i/n)
							
				// Calculate the number of jumps for each charge
				jump[event_charge[event]] += 1.0;
				
				// Set the occupancy of the grid
				grid_occ[event_mol_index[event]][curr_grid[event_charge[event]]] = false;
				
				// Set the new position in the grid from temp values
				curr_grid[event_charge[event]] = tmp_curr_grid;
				curr_mol[event_charge[event]] = tmp_curr_mol;
				
				// Check if the charge reached the end of the grid
				if ((F_dir.compare("a") == 0 && tmp_curr_grid_a >= n_mini_grid_a-1) || (F_dir.compare("b") == 0 && tmp_curr_grid_b >= n_mini_grid_b-1) || (F_dir.compare("c") == 0 && tmp_curr_grid_c >= n_mini_grid_c-1) || (F_dir.compare("ab") == 0 && fabs(dist[event_charge[event]]) > dist_tot) || (F_dir.compare("ac") == 0 && fabs(dist[event_charge[event]]) > dist_tot) || (F_dir.compare("bc") == 0 && fabs(dist[event_charge[event]]) > dist_tot)){					
					// The charge is removed
					// curr_mol.erase(curr_mol.begin()+event_charge[event]);
					// curr_grid.erase(curr_grid.begin()+event_charge[event]);
					
					// The charge appears at the beginning of the system
					int *pos;
					pos = new int[2];
					
					if (F_dir.compare("a") == 0  || F_dir.compare("b") == 0 || F_dir.compare("c") == 0){
						if (n_layer == 0)
							Dispatch_Mol_begin(i, grid_occ, pos);
						else
							Dispatch_Mol_begin_layer(i, grid_occ, pos);
					}
					else {
						if (n_layer == 0)
							Dispatch_Mol_RND(i, grid_occ, pos);
						else
							Dispatch_Mol_RND_layer(i, grid_occ, pos);
					}
					
					curr_grid[event_charge[event]] = pos[1];
					curr_mol[event_charge[event]] = pos[0];
					grid_occ[pos[0]][pos[1]] = true;
					dist[event_charge[event]] = 0.0;
					
					delete [] pos;
					
					// Print summary for the current try in the file
					if ((charge_try/n_charges) == (double(charge_try)/double(n_charges))) {
						pFile=fopen(OUT_SIMU.str().c_str(), "a");
						if (pFile==NULL) {
							int wait = 0; 
							while (wait<10 && pFile==NULL){
								cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) << " seconds to write a file" << endl;
								usleep(10*(wait+1)*(wait+1));
								pFile=fopen(OUT_SIMU.str().c_str(), "a");
								wait++;
							}
							if (wait==10 && pFile==NULL){
								cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
								exit (1);
							}
						}
						fprintf(pFile,"Frame = %d\n", i);
						fprintf(pFile,"Electric_Field_Angle = %d\n", int(F_angle));
						fprintf(pFile,"Electric_Field_Unit_Vectors: (%f, %f, %f)\n", uF_x, uF_y, uF_z);
						fprintf(pFile,"Number_of_Charges = %d\n", n_charges);
						fprintf(pFile,"Density_of_Charges = %.5e charges/cm3\n", double(n_charges)/(vol_box[i]*n_mini_grid_a*n_mini_grid_b*n_mini_grid_c*1e-24));
						fprintf(pFile,"Time_try_%d = %e\n", (charge_try/n_charges), total_time_try/(charge_try/n_charges));
						fprintf(pFile,"Distance_try_%d = %e\n", (charge_try/n_charges), total_dist_try/(charge_try/n_charges));
						fprintf(pFile,"Mu_try_%d = %lf\n", (charge_try/n_charges), total_dist_try/(total_time_try*F_norm));
						fprintf(pFile,"-------------------------------------------------------------------------------\n");
						fclose(pFile);
					}
					
					// Print info on standard output
					if ((charge_try/n_charges) % 10 == 0 && (charge_try/n_charges) == (double(charge_try)/double(n_charges))) {
						cout << "[INFO] Running Frame " << i+1 << "/" << n_frame << ", Try " << (charge_try/n_charges) << "/" << n_try << endl;
					}
					
					charge_try++;
					
				}
				
				else{
					grid_occ[tmp_curr_mol][tmp_curr_grid] = true;
				}
				
				previous_jump_ok = true;
				exit_loop = true;
				first_calc = false;
			}
		}

		// Calculates the final mobility for the frame
		mu_frame.push_back(total_dist_try/(total_time_try*F_norm));

		// Writes a summary for the frame
		pFile=fopen(OUT_SIMU.str().c_str(), "a");
		if (pFile==NULL) {
			int wait = 0; 
			while (wait<10 && pFile==NULL){
				cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) << " seconds to write a file" << endl;
				usleep(10*(wait+1)*(wait+1));
				pFile=fopen(OUT_SIMU.str().c_str(), "a");
				wait++;
			}
			if (wait==10 && pFile==NULL){
				cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
				exit (1);
			}
		}
		fprintf(pFile,"-------------------------------------------------------------------------------\n");
		fprintf(pFile,"Frame = %d\n", i);
		fprintf(pFile,"Number of tries = %d\n", n_try);
		fprintf(pFile,"Electric Field Angle = %d\n", int(F_angle));
		fprintf(pFile,"Electric Field Unit Vectors: (%f, %f, %f)\n", uF_x, uF_y, uF_z);
		fprintf(pFile,"Number of Charges = %d\n", n_charges);
		fprintf(pFile,"Density of Charges = %.5e charges/cm3\n", double(n_charges)/(vol_box[i]*n_mini_grid_a*n_mini_grid_b*n_mini_grid_c*1e-24));
		fprintf(pFile,"Average Time = %e\n", total_time_try/double(n_try));
		fprintf(pFile,"Average Distance = %e\n", total_dist_try/double(n_try));
		fprintf(pFile,"Mobility of the Frame %d = %lf\n", i, mu_frame.back());
		fclose(pFile);
		
		// ---------------------------------------------------------------------
		
		// Cleaning everything
		
		curr_mol.clear();
		curr_grid.clear();
		dist.clear(); 
		jump.clear();
		
		event_k.clear();
		event_charge.clear();
		event_mol_index.clear();
		event_neigh_num.clear();
		event_neigh_index.clear();
		
	}
	
	// Calculates the average on all the frames
	double mu_moy = 0.0; 
	double dbl_n_frame = n_frame;
	for(int i=0; i<n_frame; i++) {
		mu_moy = mu_moy + mu_frame[i];
	}
	mu_moy = mu_moy/dbl_n_frame;
	
	// Writes the final mobility
	pFile=fopen(OUT_SIMU.str().c_str(), "a");
	if (pFile==NULL) {
		int wait = 0; 
		while (wait<10 && pFile==NULL){
			cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) << " seconds to write a file" << endl;
			usleep(10*(wait+1)*(wait+1));
			pFile=fopen(OUT_SIMU.str().c_str(), "a");
			wait++;
		}
		if (wait==10 && pFile==NULL){
			cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
			exit (1);
		}
	}
	fprintf(pFile,"-------------------------------------------------\n");
	fprintf(pFile,"Number of tries = %d\n", n_try);
	fprintf(pFile,"Electric Field Angle = %d\n", int(F_angle));
	fprintf(pFile,"Electric Field Unit Vectors: (%f, %f, %f)\n", uF_x, uF_y, uF_z);
	fprintf(pFile,"Number of Charges = %d\n", n_charges);
	fprintf(pFile,"Density of Charges = %.5e charges/cm3\n", double(n_charges)/(vol_box[0]*n_mini_grid_a*n_mini_grid_b*n_mini_grid_c*1e-24));
	fprintf(pFile,"Mobility = %lf\n", mu_moy);
	fclose(pFile);
	
}

// FRM algorithm
void MC_FRM_MT(string output_folder){
	
	cout << "[INFO] Using the FRM algorithm with multithread." << endl;
	
	// Average mobility for each frame
	vector<double> mu_frame (n_frame, 0.0);
	
	double uF_x, uF_y, uF_z;
	uF_x = F_x/F_norm;
	uF_y = F_y/F_norm;
	uF_z = F_z/F_norm;

	#pragma omp parallel for private(grid_occ) //if (n_frame >= int(omp_get_max_threads()))
	
	// Start the FRM algorithm
	for (int i=0; i<n_frame; i++){
		
		// Variables for each charge
		vector<int> curr_mol, curr_grid; // Number of the molecule in the mini-grid, and number of this mini-grid
		vector<double> dist, jump; // Distance and number of jumps for each charge
		
		// Variables for a chosen event
		vector<double> event_k; // Transfer rate
		vector<int> event_charge, event_mol_index, event_neigh_num, event_neigh_index; // Initial molecule and neighbor
		
		// Total time and distance
		double total_time_try, total_dist_try;
		
		// ---------------------------------------------------------------------
		
		// Print summary for the current try
		stringstream OUT_SIMU, OUT_ERROR;
		
		OUT_SIMU << output_folder << "/simu_" << charge.c_str() << "_" << F_dir.c_str()  << "_f_" << i << ".out";
		OUT_ERROR << output_folder << "/error_" << charge.c_str() << "_" << F_dir.c_str() << "_f_" << i << ".out";
	  
		FILE * pFile;

		pFile=fopen(OUT_SIMU.str().c_str(), "w");
		if (pFile==NULL) {
			int wait = 0; 
			while (wait<10 && pFile==NULL){
				cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) << " seconds to write a file" << endl;
				usleep(10*(wait+1)*(wait+1));
				pFile=fopen(OUT_SIMU.str().c_str(), "w");
				wait++;
			}
			if (wait==10 && pFile==NULL){
				cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
				exit (1);
			}
		}
		fprintf(pFile,"===============================================================================\n");
		fprintf(pFile,"-------------------------------------------------------------------------------\n");
		fclose(pFile);	
		
		// ---------------------------------------------------------------------
		
		// Generate the grid
		grid_occ.clear();
				
		for (int ii=0; ii<n_mol; ii++){
			grid_occ.push_back( vector<bool> ());
			
			for (unsigned int x=0; x<box_neigh_label.size(); x++){
				grid_occ[ii].push_back( false );
			}
		}

		// ---------------------------------------------------------------------

		// Put the charges in the grid
		curr_mol.clear();
		curr_grid.clear();
		int *pos;
		pos = new int[2];
		
		for (unsigned int charge_i=0; charge_i<n_charges; charge_i++){
			if (n_layer == 0)
				Dispatch_Mol_RND(i, grid_occ, pos);
			else
				Dispatch_Mol_RND_layer(i, grid_occ, pos);
			
			curr_mol.push_back(pos[0]);
			curr_grid.push_back(pos[1]);
			grid_occ[pos[0]][pos[1]] = true;
		}

		delete [] pos;
		
		// ---------------------------------------------------------------------

		// Set distances, number of jumps and the travel time to zero
		total_time_try = 0.0;
		total_dist_try = 0.0;

		dist.clear();
		jump.clear(); 
		for (unsigned int charge_i=0; charge_i<n_charges; charge_i++){
			dist.push_back(0.0);
			jump.push_back(0.0);
		}

		bool previous_jump_ok = true;
		bool exit_loop = false;
		
		// ---------------------------------------------------------------------
				
		// Calculates the waiting times for all the charges
		// at the beginning of the simulation
		
		event_k.clear();
		event_charge.clear();
		event_mol_index.clear();
		event_neigh_num.clear(); 
		event_neigh_index.clear();
		
		for (unsigned int charge_i=0; charge_i<curr_mol.size(); charge_i++){
			
			double event_k_tmp = numeric_limits<double>::max();
			int event_charge_tmp = 0;
			int event_mol_index_tmp = 0;
			int event_neigh_num_tmp = 0;
			int event_neigh_index_tmp = 0;
			
			int tmp_mol_index = curr_mol[charge_i];
									
			for (unsigned int jj=0; jj<neigh_label[i][tmp_mol_index].size(); jj++){
				
				// Find index of neighbor
				int tmp_neigh_num = jj;
				int tmp_neigh_index = tmp_mol_index;
				for (int ii=0; ii<n_mol; ii++){
					if (neigh_label[i][tmp_mol_index][jj] == mol_label[ii]){
						tmp_neigh_index = ii;
						break;
					}
				}
				
				// Calculate waiting_time
				double k_tmp = numeric_limits<double>::min();
				if (tmp_neigh_index != tmp_mol_index){
					k_tmp = Marcus_Levich_Jortner_rate_electro(i, tmp_mol_index, tmp_neigh_index, tmp_neigh_num, d_x[i][tmp_mol_index][jj], d_y[i][tmp_mol_index][jj], d_z[i][tmp_mol_index][jj], dE[i][tmp_mol_index][jj], J_H[i][tmp_mol_index][jj], J_L[i][tmp_mol_index][jj], curr_mol, curr_grid, charge_i);
				
					double random_number = Rand_0_1();
					k_tmp = -log(random_number)/(k_tmp);
				
					// Keep the value only if k is larger
					if (k_tmp < event_k_tmp){
						event_k_tmp = k_tmp;
						event_charge_tmp = charge_i;
						event_mol_index_tmp = tmp_mol_index;
						event_neigh_num_tmp = jj;
						event_neigh_index_tmp = tmp_neigh_index;
					}
				}
			}
			
			// Keep only the highest k
			event_k.push_back(event_k_tmp);
			event_charge.push_back(event_charge_tmp);
			event_mol_index.push_back(event_mol_index_tmp);
			event_neigh_num.push_back(event_neigh_num_tmp);
			event_neigh_index.push_back(event_neigh_index_tmp);
		}
		
		// TO DO : SORT THE TABLE ?
		
		bool first_calc = true;
				
		// ---------------------------------------------------------------------

		unsigned int event = 0;
		int charge_previously_selected = 0;
		
		unsigned int charge_try = 0;
		while (charge_try<n_try*n_charges){	
			
			// Calculates the transfer rate for the previously selected molecule
			if (first_calc==false){

				int charge_i = charge_previously_selected;
										
				double event_k_tmp = numeric_limits<double>::max();
				int event_charge_tmp = 0;
				int event_mol_index_tmp = 0;
				int event_neigh_num_tmp = 0;
				int event_neigh_index_tmp = 0;
				
				int tmp_mol_index = curr_mol[charge_i];
										
				for (unsigned int jj=0; jj<neigh_label[i][tmp_mol_index].size(); jj++){
					
					// Find index of neighbor
					int tmp_neigh_num = jj;
					int tmp_neigh_index = tmp_mol_index;
					for (int ii=0; ii<n_mol; ii++){
						if (neigh_label[i][tmp_mol_index][jj] == mol_label[ii]){
							tmp_neigh_index = ii;
							break;
						}
					}
					
					// Calculate waiting_time
					double k_tmp = numeric_limits<double>::max();
					if (tmp_neigh_index != tmp_mol_index){
						k_tmp = Marcus_Levich_Jortner_rate_electro(i, tmp_mol_index, tmp_neigh_index, tmp_neigh_num, d_x[i][tmp_mol_index][jj], d_y[i][tmp_mol_index][jj], d_z[i][tmp_mol_index][jj], dE[i][tmp_mol_index][jj], J_H[i][tmp_mol_index][jj], J_L[i][tmp_mol_index][jj], curr_mol, curr_grid, charge_i);
					
						double random_number = Rand_0_1();
						k_tmp = -log(random_number)/(k_tmp);
				
						// Keep the value only if k is larger
						if (k_tmp < event_k_tmp){
							event_k_tmp = k_tmp;
							event_charge_tmp = charge_i;
							event_mol_index_tmp = tmp_mol_index;
							event_neigh_num_tmp = jj;
							event_neigh_index_tmp = tmp_neigh_index;
						}
					}
				}
				// Keep the highest k
				event_k[charge_previously_selected] = event_k_tmp;
				event_charge[charge_previously_selected] = event_charge_tmp;
				event_mol_index[charge_previously_selected] = event_mol_index_tmp;
				event_neigh_num[charge_previously_selected] = event_neigh_num_tmp;
				event_neigh_index[charge_previously_selected] = event_neigh_index_tmp;
			}
			
			else{
				exit_loop = false;
			}
			
			// Choose the fastest event
			double k_tmp = numeric_limits<double>::max();
			for (unsigned int t=0; t<event_k.size(); t++){
				//cout << event_k[t] << endl;
				if(event_k[t] < k_tmp) {
					k_tmp = event_k[t];
					charge_previously_selected = event_charge[t];
					event = t;
				}
			}

			// =======================================
			//
			// Summary of variables:
			// ---------------------
			//
			// event_charge[event]: index of the charge of the event (= charge_i)
			// event_mol_index[event]: index of the molecule where the chosen charge is (= ii)
			// event_neigh_num[event]: number of the neighbor of the molecule where the chosen charge is (= jj)
			// event_neigh_index[event]: index of the molecule where the chosen charge is going to jump (= ii for the neighbor)
			// curr_mol[event_charge[event]]: index of the molecule where the chosen charge is (= ii)
			// curr_grid[event_charge[event]]: index of the mini-grid where the chosen charge is (= x)
			//
			// =======================================
	
			// Calculate the new position in the grid
			int tmp_curr_mol = 0; // Expected final molecule index
			int tmp_curr_grid = 0; // Expected final mini-grid index
			int tmp_curr_grid_a = 0, tmp_curr_grid_b = 0, tmp_curr_grid_c = 0;
			bool out_of_system = true;
			
			// Find next mini-grid
			if (neigh_jump_vec_a[i][event_mol_index[event]][event_neigh_num[event]] == 0 && neigh_jump_vec_b[i][event_mol_index[event]][event_neigh_num[event]] == 0 && neigh_jump_vec_c[i][event_mol_index[event]][event_neigh_num[event]] == 0){
				tmp_curr_grid = curr_grid[event_charge[event]];
				tmp_curr_grid_a = box_a[tmp_curr_grid];
				tmp_curr_grid_b = box_b[tmp_curr_grid];
				tmp_curr_grid_c = box_c[tmp_curr_grid];
				out_of_system = false;
			}
			
			else {
				for (unsigned int xx=0; xx<box_neigh_label[curr_grid[event_charge[event]]].size(); xx++){
					if (neigh_jump_vec_a[i][event_mol_index[event]][event_neigh_num[event]] == box_neigh_a[curr_grid[event_charge[event]]][xx] && neigh_jump_vec_b[i][event_mol_index[event]][event_neigh_num[event]] == box_neigh_b[curr_grid[event_charge[event]]][xx] && neigh_jump_vec_c[i][event_mol_index[event]][event_neigh_num[event]] == box_neigh_c[curr_grid[event_charge[event]]][xx]){
						tmp_curr_grid = box_neigh_label[curr_grid[event_charge[event]]][xx];
						tmp_curr_grid_a = box_a[tmp_curr_grid];
						tmp_curr_grid_b = box_b[tmp_curr_grid];
						tmp_curr_grid_c = box_c[tmp_curr_grid];
						out_of_system = false;
						break;
					}
				}
			}
				
			// Find next molecule
			tmp_curr_mol = event_neigh_index[event];
			
			// Do not jump if the charge goes out of the system
			if (out_of_system){
				previous_jump_ok = false;
				exit_loop = true;
				first_calc = false;
				//printf("out of system\n");
			}
			
			// Do not jump if the next molecule is occupied
			else if (grid_occ[tmp_curr_mol][tmp_curr_grid] == true){
				previous_jump_ok = false;
				exit_loop = true;
				first_calc = false;
				//printf("occupied\n");
			}

			// Do not jump if the charge change of layer
			else if (n_layer > 0 && list_layer[i][curr_mol[event_charge[event]]] != list_layer[i][tmp_curr_mol]){
				previous_jump_ok = false;
				exit_loop = true;
				first_calc = false;
			}
			
			else{
				
				// Calculate the total time
				total_time_try += event_k[event];
	
				// Substract the waiting time
				double k_event = event_k[event];
				for (unsigned int t=0; t<event_k.size(); t++){
					event_k[t] = event_k[t] - k_event;
				}

				// Calculate the distance traveled by the charge and the total distance
				double event_dist = (d_x[i][event_mol_index[event]][event_neigh_num[event]] * uF_x + d_y[i][event_mol_index[event]][event_neigh_num[event]] * uF_y + d_z[i][event_mol_index[event]][event_neigh_num[event]] * uF_z)*1E-8;
				dist[event_charge[event]] += event_dist;
				//total_dist_try += event_dist;
				total_dist_try += event_dist/double(curr_mol.size());
				// Note:
				// We divide here by the number of charges in the system, 
				// so there is no problem afterwards if the number of charges
				// in the system change. Indeed:
				// <x> = (1/n) * sum_i x_i = sum_i (x_i/n)
							
				// Calculate the number of jumps for each charge
				jump[event_charge[event]] += 1.0;
				
				// Set the occupancy of the grid
				grid_occ[event_mol_index[event]][curr_grid[event_charge[event]]] = false;
				
				// Set the new position in the grid from temp values
				curr_grid[event_charge[event]] = tmp_curr_grid;
				curr_mol[event_charge[event]] = tmp_curr_mol;

				// Check if the charge reached the end of the grid
				if ((F_dir.compare("a") == 0 && tmp_curr_grid_a >= n_mini_grid_a-1) || (F_dir.compare("b") == 0 && tmp_curr_grid_b >= n_mini_grid_b-1) || (F_dir.compare("c") == 0 && tmp_curr_grid_c >= n_mini_grid_c-1) || (F_dir.compare("ab") == 0 && fabs(dist[event_charge[event]]) > dist_tot) || (F_dir.compare("ac") == 0 && fabs(dist[event_charge[event]]) > dist_tot) || (F_dir.compare("bc") == 0 && fabs(dist[event_charge[event]]) > dist_tot)){
				
					// The charge is removed
					// curr_mol.erase(curr_mol.begin()+event_charge[event]);
					// curr_grid.erase(curr_grid.begin()+event_charge[event]);
					
					// The charge appears at the beginning of the system
					int *pos;
					pos = new int[2];
					
					if (F_dir.compare("a") == 0  || F_dir.compare("b") == 0 || F_dir.compare("c") == 0){
						if (n_layer == 0)
							Dispatch_Mol_begin(i, grid_occ, pos);
						else
							Dispatch_Mol_begin_layer(i, grid_occ, pos);
					}
					else {
						if (n_layer == 0)
							Dispatch_Mol_RND(i, grid_occ, pos);
						else
							Dispatch_Mol_RND_layer(i, grid_occ, pos);
					}
					
					curr_grid[event_charge[event]] = pos[1];
					curr_mol[event_charge[event]] = pos[0];
					grid_occ[pos[0]][pos[1]] = true;
					dist[event_charge[event]] = 0.0;
					
					delete [] pos;
					
					// Print summary for the current try in the file
					if ((charge_try/n_charges) == (double(charge_try)/double(n_charges))) {
						pFile=fopen(OUT_SIMU.str().c_str(), "a");
						if (pFile==NULL) {
							int wait = 0; 
							while (wait<10 && pFile==NULL){
								cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) << " seconds to write a file" << endl;
								usleep(10*(wait+1)*(wait+1));
								pFile=fopen(OUT_SIMU.str().c_str(), "a");
								wait++;
							}
							if (wait==10 && pFile==NULL){
								cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
								exit (1);
							}
						}
						fprintf(pFile,"Frame = %d\n", i);
						fprintf(pFile,"Electric_Field_Angle = %d\n", int(F_angle));
						fprintf(pFile,"Electric_Field_Unit_Vectors: (%f, %f, %f)\n", uF_x, uF_y, uF_z);
						fprintf(pFile,"Number_of_Charges = %d\n", n_charges);
						fprintf(pFile,"Density_of_Charges = %.5e charges/cm3\n", double(n_charges)/(vol_box[i]*n_mini_grid_a*n_mini_grid_b*n_mini_grid_c*1e-24));
						fprintf(pFile,"Time_try_%d = %e\n", (charge_try/n_charges), total_time_try/(charge_try/n_charges));
						fprintf(pFile,"Distance_try_%d = %e\n", (charge_try/n_charges), total_dist_try/(charge_try/n_charges));
						fprintf(pFile,"Mu_try_%d = %lf\n", (charge_try/n_charges), total_dist_try/(total_time_try*F_norm));
						fprintf(pFile,"-------------------------------------------------------------------------------\n");
						fclose(pFile);
					}
					
					// Print info on standard output
					if ((charge_try/n_charges) % 10 == 0 && (charge_try/n_charges) == (double(charge_try)/double(n_charges))) {
						cout << "[INFO] Running Frame " << i+1 << "/" << n_frame << ", Try " << (charge_try/n_charges) << "/" << n_try << endl;
					}
					
					charge_try++;
					
				}
				
				else{
					grid_occ[tmp_curr_mol][tmp_curr_grid] = true;
				}
				
				previous_jump_ok = true;
				exit_loop = true;
				first_calc = false;
			}
		}

		// Calculates the final mobility for the frame
		mu_frame[i] = total_dist_try/(total_time_try*F_norm);

		// Writes a summary for the frame
		pFile=fopen(OUT_SIMU.str().c_str(), "a");
		if (pFile==NULL) {
			int wait = 0; 
			while (wait<10 && pFile==NULL){
				cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) << " seconds to write a file" << endl;
				usleep(10*(wait+1)*(wait+1));
				pFile=fopen(OUT_SIMU.str().c_str(), "a");
				wait++;
			}
			if (wait==10 && pFile==NULL){
				cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
				exit (1);
			}
		}
		fprintf(pFile,"-------------------------------------------------------------------------------\n");
		fprintf(pFile,"Frame = %d\n", i);
		fprintf(pFile,"Number of tries = %d\n", n_try);
		fprintf(pFile,"Electric Field Angle = %d\n", int(F_angle));
		fprintf(pFile,"Electric Field Unit Vectors: (%f, %f, %f)\n", uF_x, uF_y, uF_z);
		fprintf(pFile,"Number of Charges = %d\n", n_charges);
		fprintf(pFile,"Density of Charges = %.5e charges/cm3\n", double(n_charges)/(vol_box[i]*n_mini_grid_a*n_mini_grid_b*n_mini_grid_c*1e-24));
		fprintf(pFile,"Average Time = %e\n", total_time_try/double(n_try));
		fprintf(pFile,"Average Distance = %e\n", total_dist_try/double(n_try));
		fprintf(pFile,"Mobility = %lf\n", mu_frame[i]);
		fclose(pFile);
		
		// ---------------------------------------------------------------------
		
		// Cleaning everything
		
		curr_mol.clear();
		curr_grid.clear();
		dist.clear(); 
		jump.clear();
		
		event_k.clear();
		event_charge.clear();
		event_mol_index.clear();
		event_neigh_num.clear();
		event_neigh_index.clear();
		
	}
	
	#pragma omp barrier
	
	// Calculates the average on all the frames
	double mu_moy = 0.0; 
	double dbl_n_frame = n_frame;
	for(int i=0; i<n_frame; i++) {
		mu_moy = mu_moy + mu_frame[i];
	}
	mu_moy = mu_moy/dbl_n_frame;
	
	// Writes the final mobility
	stringstream OUT_SIMU, OUT_ERROR;
	
	OUT_SIMU << output_folder << "/simu_" << charge.c_str() << "_" << F_dir.c_str() << ".out";
	OUT_ERROR << output_folder << "/error_" << charge.c_str() << "_" << F_dir.c_str() << ".out";
	  
	FILE * pFile;
		
	pFile=fopen(OUT_SIMU.str().c_str(), "a");
	if (pFile==NULL) {
		int wait = 0; 
		while (wait<10 && pFile==NULL){
			cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) << " seconds to write a file" << endl;
			usleep(10*(wait+1)*(wait+1));
			pFile=fopen(OUT_SIMU.str().c_str(), "a");
			wait++;
		}
		if (wait==10 && pFile==NULL){
			cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
			exit (1);
		}
	}
	fprintf(pFile,"-------------------------------------------------\n");
	fprintf(pFile,"Number of tries = %d\n", n_try);
	fprintf(pFile,"Electric Field Angle = %d\n", int(F_angle));
	fprintf(pFile,"Electric Field Unit Vectors: (%f, %f, %f)\n", uF_x, uF_y, uF_z);
	fprintf(pFile,"Number of Charges = %d\n", n_charges);
	fprintf(pFile,"Density of Charges = %.5e charges/cm3\n", double(n_charges)/(vol_box[0]*n_mini_grid_a*n_mini_grid_b*n_mini_grid_c*1e-24));
	fprintf(pFile,"Mobility = %lf\n", mu_moy);
	fclose(pFile);
	
}

int main(int argc, char **argv){
	
	cout << "\n===============================================================================" << endl;
	cout << "-------------------------- Starting a KMC simulation --------------------------" << endl;
	cout << "===============================================================================" << endl << endl;
	
	time_t t_start, t_stop;
	string timeinfo_start, timeinfo_stop;
	
	// Get start time
	t_start = time(NULL);
	timeinfo_start = asctime(localtime(&t_start));
	
	srand(time(NULL));

	// Define specific varables for layer case
	n_layer = 4;
	layer = numeric_limits<int>::max();
	min_layer.assign(n_layer, 0.0);
	max_layer.assign(n_layer, 0.0);

	int s;
	string input_file;
	string input_folder = ".";
	string output_folder = ".";
	string method;
	
  	while ((s = getopt_long (argc, argv, "I:i:o:c:d:n:m:l:", NULL, NULL)) != -1){
      	switch (s){
			case 'I':
				input_file = optarg;
				break;
				
			case 'i':
				input_folder = optarg;
				struct stat st_in;
				if(stat(input_folder.c_str(), &st_in) != 0){
					cout << "[ERROR] Input folder " << input_folder << " doesn't exist! Exiting..." << endl;
					Clear_All();
					exit(1);
				}
	  			break;
	  			
			case 'o':
				output_folder = optarg;
				struct stat st_out;
				if(stat(output_folder.c_str(), &st_out) != 0){
					mkdir(output_folder.c_str(), 0750);
					cout << "[INFO] Output folder " << output_folder << " created." << endl;
				}
	  			break;
	  			
			case 'c':
				charge = optarg;
	  			break;
	  			
	  		case 'd':
				F_dir = optarg;
				break;
	  			
	  		case 'n':
				n_charges = atoi(optarg);
				break;
	  			
	  		case 'm':
				method = optarg;
	  			break;
	  		
			case 'l':
				layer = atoi(optarg);
				break;
			}
	}
	
	// Check that the charge and the direction of the electric field are specified
	if (charge.compare("e") == 0 || charge.compare("h") == 0)
		cout << "[INFO] The charge is: " << charge << endl;
		
	else {
		cerr << "[ERROR] Charge not specified! Please use -c {e,h}. Exiting..." << endl;
		Clear_All();
		exit(1);
	}
	
	if (F_dir.compare("a") == 0 || F_dir.compare("b") == 0 || F_dir.compare("c") == 0)
		cout << "[INFO] Electric field is along the '" << F_dir << "' direction." << endl;
		
	else if (F_dir.compare("ab") == 0 || F_dir.compare("ac") == 0 || F_dir.compare("bc") == 0)
		cout << "[INFO] Electric field will probe the anisotropy in the '" << F_dir << "'plane." << endl;
		
	else {
		cerr << "[ERROR] Direction of the electric field not specified! Please use -d {a,b,c,ab,ac,bc}. Exiting..." << endl;
		Clear_All();
		exit(1);
	}

	if (layer != std::numeric_limits<int>::max())
		cout << "[INFO] The charges will be set in layer " << layer << endl;
		
	else {
		cerr << "[ERROR] Layer not specified. Please use -l. Exiting..." << endl;
		Clear_All();
		exit(1);
	}

	// Read the required files
	Read_MC(input_file, input_folder, false);
	Read_CELL(input_file, input_folder, false);
	Read_CM(input_file, input_folder, false);
	Read_E_av(input_file, input_folder, false);

	// Set up lambda_i
	if (charge.compare("e") == 0) {
		LAMBDA_I = LAMBDA_I_E;
	}
	else {
		LAMBDA_I = LAMBDA_I_H;
	}

	// Calculates Electric field vectors
	Calcul_F_vector(false);
	
	// Calculates distances and DeltaE
	Calcul_Dist(false);
	Calcul_DeltaE(false);

	// Build the grid
	Find_Layer(false);
	Build_Grid(false);
	
	// Save the triangular matrix
	vector< vector< vector<int> > > neigh_label_ref = neigh_label;
	vector< vector< vector<int> > > neigh_jump_vec_a_ref = neigh_jump_vec_a, neigh_jump_vec_b_ref = neigh_jump_vec_b, neigh_jump_vec_c_ref = neigh_jump_vec_c;
	vector< vector< vector<double> > > J_H_ref = J_H, J_L_ref = J_L;
	vector< vector< vector<double> > > d_x_ref = d_x, d_y_ref = d_y, d_z_ref = d_z;
	vector< vector< vector<double> > > dE_ref = dE;				
	
	//#pragma omp parallel for private(F_x, F_y, F_z, F_angle, F_x_list, F_y_list, F_z_list, F_angle_list, k, k_inv, grid_occ) //if (n_frame < int(omp_get_max_threads()) && anisotropy == true)
	for (unsigned int m=0; m<F_angle_list.size(); m++) {
		
		F_x = F_x_list[m];
		F_y = F_y_list[m];
		F_z = F_z_list[m];
		F_angle = F_angle_list[m];

		// Calculate transfer rates and the full matrix, mostly for information
		Marcus_Levich_Jortner_CST(); // Calculates constants
		Calcul_k(false);
		Full_Matrix();

		// Print a summary
		Print_Summary(output_folder);
		
		// Calculates 1/k and clear the k table
		Inverse_Clear_k(false);
		k.clear();

		if(method.compare("bkl") == 0) {
			if(MT) {
				MC_BKL_MT(output_folder);
			}
			else {
				MC_BKL(output_folder);
			}
		}

		else {
			if(MT) {
				MC_FRM_MT(output_folder);
			}
			else {
				MC_FRM(output_folder);
			}
		}

		// Clear everything (not necessary) and set to reference values
		neigh_label.clear(); neigh_label = neigh_label_ref;
		neigh_jump_vec_a.clear(); neigh_jump_vec_a = neigh_jump_vec_a_ref;
		neigh_jump_vec_b.clear(); neigh_jump_vec_b = neigh_jump_vec_b_ref;
		neigh_jump_vec_c.clear(); neigh_jump_vec_c = neigh_jump_vec_c_ref;
		J_H.clear(); J_H = J_H_ref;
		J_L.clear(); J_L = J_L_ref;
		d_x.clear(); d_x = d_x_ref;
		d_y.clear(); d_y = d_y_ref;
		d_z.clear(); d_z = d_z_ref;
		dE.clear();	dE = dE_ref;
		k.clear();
		k_inv.clear();
	}
	//#pragma omp barrier

	Clear_All();
		
	// Get stop time
	t_stop = time(NULL);
	timeinfo_stop = asctime(localtime(&t_stop));

	// Print final information
	cout << "[INFO] Start time: " << timeinfo_start;	
	cout << "[INFO] Stop time: " << timeinfo_stop;
	cout << "[INFO] Calculation took " << t_stop-t_start << " seconds." << endl << endl;
	cout << "===============================================================================" << endl;
	cout << "-------------------- Exiting normally, everything was ok! ---------------------" << endl;
	cout << "===============================================================================\n" << endl;
	
	return 0;
}
