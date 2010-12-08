/**
 *******************************************************************************
 * Copyright (C) 2010 Nicolas Martinelli, nicolas.martinelli@gmail.com         *
 * This library is part of the MC_BKL and MC_BKL_layer software                *
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

// C++ libraries
#include <iostream>
#include <fstream> 
#include <string> 
#include <limits> 
#include <sstream> 
#include <vector>

// C libraries
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include <omp.h>

// Local libraries
#include "variables.h"
#include "clear.h"
#include "mathplus.h"
#include "transferrate.h"
#include "dispatch.h"
#include "printsummary.h"

using namespace std;

// FRM algorithm
void MC_FRM(string output_folder){
	
	cout << "[INFO] Using the FRM algorithm." << endl;
	
	// Variables for each charge
	vector<int> curr_mol, curr_box; // Number of the molecule in the mini-grid, and number of this mini-grid
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
	
	// Start the FRM algorithm
	for (int i=0; i<n_frame; i++){
		
		// ---------------------------------------------------------------------
		
		// Set up variables for output files
		stringstream OUT_SIMU_FRAME, OUT_ERROR;
		
		OUT_SIMU_FRAME << output_folder << "/simu_" << charge.c_str() << "_" << F_dir.c_str() << ".out";
		OUT_ERROR << output_folder << "/error_" << charge.c_str() << "_" << F_dir.c_str() << ".out";

		// ---------------------------------------------------------------------
		
		// Generate the grid
		grid_occ.clear();
				
		for (unsigned int x=0; x<box_neigh_label.size(); x++){
			grid_occ.push_back( vector<bool> ());
			
			for (int ii=0; ii<n_mol; ii++){
				grid_occ[x].push_back( false );
			}
		}

		// ---------------------------------------------------------------------

		// Put the charges in the grid
		curr_mol.clear();
		curr_box.clear();
		int *pos;
		pos = new int[2];
		
		for (unsigned int charge_i=0; charge_i<n_charges; charge_i++){
			if (n_layer == 0)
				Dispatch_Mol_RND(i, grid_occ, pos);
			else
				Dispatch_Mol_RND_layer(i, grid_occ, pos);
			
			curr_mol.push_back(pos[0]);
			curr_box.push_back(pos[1]);
			grid_occ[pos[1]][pos[0]] = true;
			//grid_probability[i][pos[1]][pos[0]][charge_i] += 1.0;
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
					k_tmp = Marcus_Levich_Jortner_rate_electro(i, tmp_mol_index, tmp_neigh_index, tmp_neigh_num, d_x[i][tmp_mol_index][jj], d_y[i][tmp_mol_index][jj], d_z[i][tmp_mol_index][jj], dE[i][tmp_mol_index][jj], J_H[i][tmp_mol_index][jj], J_L[i][tmp_mol_index][jj], curr_mol, curr_box, charge_i);
				
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
						k_tmp = Marcus_Levich_Jortner_rate_electro(i, tmp_mol_index, tmp_neigh_index, tmp_neigh_num, d_x[i][tmp_mol_index][jj], d_y[i][tmp_mol_index][jj], d_z[i][tmp_mol_index][jj], dE[i][tmp_mol_index][jj], J_H[i][tmp_mol_index][jj], J_L[i][tmp_mol_index][jj], curr_mol, curr_box, charge_i);
					
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
			// curr_box[event_charge[event]]: index of the mini-grid where the chosen charge is (= x)
			//
			// =======================================
	
			// Calculate the new position in the grid
			int tmp_curr_mol = 0; // Expected final molecule index
			int tmp_curr_box = 0; // Expected final mini-grid index
			int tmp_curr_box_a = 0, tmp_curr_box_b = 0, tmp_curr_box_c = 0;
			bool out_of_system = true;
			
			// Find next mini-grid
			if (neigh_jump_vec_a[i][event_mol_index[event]][event_neigh_num[event]] == 0 && neigh_jump_vec_b[i][event_mol_index[event]][event_neigh_num[event]] == 0 && neigh_jump_vec_c[i][event_mol_index[event]][event_neigh_num[event]] == 0){
				tmp_curr_box = curr_box[event_charge[event]];
				tmp_curr_box_a = box_a[tmp_curr_box];
				tmp_curr_box_b = box_b[tmp_curr_box];
				tmp_curr_box_c = box_c[tmp_curr_box];
				out_of_system = false;
				
			}
			
			else {
				for (unsigned int xx=0; xx<box_neigh_label[curr_box[event_charge[event]]].size(); xx++){
					if (neigh_jump_vec_a[i][event_mol_index[event]][event_neigh_num[event]] == box_neigh_a[curr_box[event_charge[event]]][xx] && neigh_jump_vec_b[i][event_mol_index[event]][event_neigh_num[event]] == box_neigh_b[curr_box[event_charge[event]]][xx] && neigh_jump_vec_c[i][event_mol_index[event]][event_neigh_num[event]] == box_neigh_c[curr_box[event_charge[event]]][xx]){
						tmp_curr_box = box_neigh_label[curr_box[event_charge[event]]][xx];
						tmp_curr_box_a = box_a[tmp_curr_box];
						tmp_curr_box_b = box_b[tmp_curr_box];
						tmp_curr_box_c = box_c[tmp_curr_box];
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
			else if (grid_occ[tmp_curr_box][tmp_curr_mol] == true){
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
				grid_occ[curr_box[event_charge[event]]][event_mol_index[event]] = false;
				grid_probability[i][curr_box[event_charge[event]]][event_mol_index[event]][event_charge[event]] += event_dist;
				
				// Set the new position in the grid from temp values
				curr_box[event_charge[event]] = tmp_curr_box;
				curr_mol[event_charge[event]] = tmp_curr_mol;
				
				// Check if the charge reached the end of the grid
				if ((F_dir.compare("a") == 0 && tmp_curr_box_a >= n_box_a-1) || (F_dir.compare("b") == 0 && tmp_curr_box_b >= n_box_b-1) || (F_dir.compare("c") == 0 && tmp_curr_box_c >= n_box_c-1) || (F_dir.compare("ab") == 0 && fabs(dist[event_charge[event]]) > dist_tot) || (F_dir.compare("ac") == 0 && fabs(dist[event_charge[event]]) > dist_tot) || (F_dir.compare("bc") == 0 && fabs(dist[event_charge[event]]) > dist_tot)){					
					// The charge is removed
					// curr_mol.erase(curr_mol.begin()+event_charge[event]);
					// curr_box.erase(curr_box.begin()+event_charge[event]);
					
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
					
					curr_box[event_charge[event]] = pos[1];
					curr_mol[event_charge[event]] = pos[0];
					grid_occ[pos[1]][pos[0]] = true;
					dist[event_charge[event]] = 0.0;
					
					delete [] pos;
					
					// Print summary for the current try in the file
					if ((charge_try/n_charges) == (double(charge_try)/double(n_charges))) {
						Print_Summary_Try(output_folder, i, charge_try, total_dist_try, total_time_try);
					}
					
					// Print info on standard output
					if ((charge_try/n_charges) % 10 == 0 && (charge_try/n_charges) == (double(charge_try)/double(n_charges))) {
						cout << "[INFO] Running Frame " << i+1 << "/" << n_frame << ", Try " << (charge_try/n_charges) << "/" << n_try << endl;
					}
					
					charge_try++;
					
				}
				
				else{
					grid_occ[tmp_curr_box][tmp_curr_mol] = true;
				}
				
				previous_jump_ok = true;
				exit_loop = true;
				first_calc = false;
			}
		}

		// Calculates the final mobility for the frame
		mu_frame.push_back(total_dist_try/(total_time_try*F_norm));

		// Writes a summary for the frame
		Print_Summary_Frame(output_folder, i, total_dist_try, total_time_try, mu_frame);
		
		// ---------------------------------------------------------------------
		
		// Cleaning everything
		
		curr_mol.clear();
		curr_box.clear();
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
	Print_Summary_Final(output_folder, mu_moy);
	
}

// FRM algorithm with multithread support
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
		vector<int> curr_mol, curr_box; // Number of the molecule in the mini-grid, and number of this mini-grid
		vector<double> dist, jump; // Distance and number of jumps for each charge
		
		// Variables for a chosen event
		vector<double> event_k; // Transfer rate
		vector<int> event_charge, event_mol_index, event_neigh_num, event_neigh_index; // Initial molecule and neighbor
		
		// Total time and distance
		double total_time_try, total_dist_try;
		
		// ---------------------------------------------------------------------
		
		// Set up variables for output files
		stringstream OUT_SIMU_FRAME, OUT_ERROR;
		
		OUT_SIMU_FRAME << output_folder << "/simu_" << charge.c_str() << "_" << F_dir.c_str()  << "_f_" << i << ".out";
		OUT_ERROR << output_folder << "/error_" << charge.c_str() << "_" << F_dir.c_str() << "_f_" << i << ".out";
		
		// ---------------------------------------------------------------------
		
		// Generate the grid
		grid_occ.clear();
				
		for (unsigned int x=0; x<box_neigh_label.size(); x++){
			grid_occ.push_back( vector<bool> ());
			
			for (int ii=0; ii<n_mol; ii++){
				grid_occ[x].push_back( false );
			}
		}

		// ---------------------------------------------------------------------

		// Put the charges in the grid
		curr_mol.clear();
		curr_box.clear();
		int *pos;
		pos = new int[2];
		
		for (unsigned int charge_i=0; charge_i<n_charges; charge_i++){
			if (n_layer == 0)
				Dispatch_Mol_RND(i, grid_occ, pos);
			else
				Dispatch_Mol_RND_layer(i, grid_occ, pos);
			
			curr_mol.push_back(pos[0]);
			curr_box.push_back(pos[1]);
			grid_occ[pos[1]][pos[0]] = true;
			//grid_probability[i][pos[1]][pos[0]][charge_i] += 1.0;
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
					k_tmp = Marcus_Levich_Jortner_rate_electro(i, tmp_mol_index, tmp_neigh_index, tmp_neigh_num, d_x[i][tmp_mol_index][jj], d_y[i][tmp_mol_index][jj], d_z[i][tmp_mol_index][jj], dE[i][tmp_mol_index][jj], J_H[i][tmp_mol_index][jj], J_L[i][tmp_mol_index][jj], curr_mol, curr_box, charge_i);
				
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
						k_tmp = Marcus_Levich_Jortner_rate_electro(i, tmp_mol_index, tmp_neigh_index, tmp_neigh_num, d_x[i][tmp_mol_index][jj], d_y[i][tmp_mol_index][jj], d_z[i][tmp_mol_index][jj], dE[i][tmp_mol_index][jj], J_H[i][tmp_mol_index][jj], J_L[i][tmp_mol_index][jj], curr_mol, curr_box, charge_i);
					
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
			// curr_box[event_charge[event]]: index of the mini-grid where the chosen charge is (= x)
			//
			// =======================================
	
			// Calculate the new position in the grid
			int tmp_curr_mol = 0; // Expected final molecule index
			int tmp_curr_box = 0; // Expected final mini-grid index
			int tmp_curr_box_a = 0, tmp_curr_box_b = 0, tmp_curr_box_c = 0;
			bool out_of_system = true;
			
			// Find next mini-grid
			if (neigh_jump_vec_a[i][event_mol_index[event]][event_neigh_num[event]] == 0 && neigh_jump_vec_b[i][event_mol_index[event]][event_neigh_num[event]] == 0 && neigh_jump_vec_c[i][event_mol_index[event]][event_neigh_num[event]] == 0){
				tmp_curr_box = curr_box[event_charge[event]];
				tmp_curr_box_a = box_a[tmp_curr_box];
				tmp_curr_box_b = box_b[tmp_curr_box];
				tmp_curr_box_c = box_c[tmp_curr_box];
				out_of_system = false;
			}
			
			else {
				for (unsigned int xx=0; xx<box_neigh_label[curr_box[event_charge[event]]].size(); xx++){
					if (neigh_jump_vec_a[i][event_mol_index[event]][event_neigh_num[event]] == box_neigh_a[curr_box[event_charge[event]]][xx] && neigh_jump_vec_b[i][event_mol_index[event]][event_neigh_num[event]] == box_neigh_b[curr_box[event_charge[event]]][xx] && neigh_jump_vec_c[i][event_mol_index[event]][event_neigh_num[event]] == box_neigh_c[curr_box[event_charge[event]]][xx]){
						tmp_curr_box = box_neigh_label[curr_box[event_charge[event]]][xx];
						tmp_curr_box_a = box_a[tmp_curr_box];
						tmp_curr_box_b = box_b[tmp_curr_box];
						tmp_curr_box_c = box_c[tmp_curr_box];
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
			else if (grid_occ[tmp_curr_box][tmp_curr_mol] == true){
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
				grid_occ[curr_box[event_charge[event]]][event_mol_index[event]] = false;
				grid_probability[i][curr_box[event_charge[event]]][event_mol_index[event]][event_charge[event]] += event_dist;
				
				// Set the new position in the grid from temp values
				curr_box[event_charge[event]] = tmp_curr_box;
				curr_mol[event_charge[event]] = tmp_curr_mol;

				// Check if the charge reached the end of the grid
				if ((F_dir.compare("a") == 0 && tmp_curr_box_a >= n_box_a-1) || (F_dir.compare("b") == 0 && tmp_curr_box_b >= n_box_b-1) || (F_dir.compare("c") == 0 && tmp_curr_box_c >= n_box_c-1) || (F_dir.compare("ab") == 0 && fabs(dist[event_charge[event]]) > dist_tot) || (F_dir.compare("ac") == 0 && fabs(dist[event_charge[event]]) > dist_tot) || (F_dir.compare("bc") == 0 && fabs(dist[event_charge[event]]) > dist_tot)){
				
					// The charge is removed
					// curr_mol.erase(curr_mol.begin()+event_charge[event]);
					// curr_box.erase(curr_box.begin()+event_charge[event]);
					
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
					
					curr_box[event_charge[event]] = pos[1];
					curr_mol[event_charge[event]] = pos[0];
					grid_occ[pos[1]][pos[0]] = true;
					dist[event_charge[event]] = 0.0;
					
					delete [] pos;
					
					// Print summary for the current try in the file
					if ((charge_try/n_charges) == (double(charge_try)/double(n_charges))) {
						Print_Summary_Try(output_folder, i, charge_try, total_dist_try, total_time_try);
					}
					
					// Print info on standard output
					if ((charge_try/n_charges) % 10 == 0 && (charge_try/n_charges) == (double(charge_try)/double(n_charges))) {
						cout << "[INFO] Running Frame " << i+1 << "/" << n_frame << ", Try " << (charge_try/n_charges) << "/" << n_try << endl;
					}
					
					charge_try++;
					
				}
				
				else{
					grid_occ[tmp_curr_box][tmp_curr_mol] = true;
				}
				
				previous_jump_ok = true;
				exit_loop = true;
				first_calc = false;
			}
		}

		// Calculates the final mobility for the frame
		mu_frame[i] = total_dist_try/(total_time_try*F_norm);

		// Writes a summary for the frame
		Print_Summary_Frame(output_folder, i, total_dist_try, total_time_try, mu_frame);
		
		// ---------------------------------------------------------------------
		
		// Cleaning everything
		
		curr_mol.clear();
		curr_box.clear();
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
	Print_Summary_Final(output_folder, mu_moy);
	
}
