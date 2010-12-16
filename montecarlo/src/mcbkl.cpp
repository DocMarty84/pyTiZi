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
#include <ctime>
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

#include <boost/random.hpp>
#include <omp.h>

// Local libraries
#include "variables.h"
#include "clear.h"
#include "mathplus.h"
#include "transferrate.h"
#include "dispatch.h"
#include "printsummary.h"
#include "deltav.h"

using namespace std;

// BKL algorithm with multithread support
void MC_BKL_MT(string output_folder){
	
	t_info = time(NULL); t_info_str = asctime(localtime(&t_info)); 
	t_info_str.erase(t_info_str.length()-1,1); 
	cout << "[INFO: " << t_info_str << "] Using the BKL algorithm with multithreading." << endl;
	
	// Create a Mersenne twister random number generator
	// that is seeded once with #seconds since 1970
	static boost::mt19937 rng(static_cast<unsigned> (time(0)));
	
	// Average mobility for each frame
	vector<double> mu_frame (n_frame, 0.0);
	
	// Electric Field
	double uF_x, uF_y, uF_z;
	uF_x = F_x/F_norm;
	uF_y = F_y/F_norm;
	uF_z = F_z/F_norm;
	
	// Charge properties
	chrg_E_electrostatic.clear();
	chrg_E_0.clear();
	chrg_E_1.clear();
	chrg_total_time.clear();
	chrg_total_dist.clear();
	
	for (int i=0; i<n_frame; i++){
		chrg_E_electrostatic.push_back( vector< vector< double > > () ); 
		chrg_E_0.push_back( vector< vector< double > > () );
		chrg_E_1.push_back( vector< vector< double > > () );
		chrg_total_time.push_back( vector< double > () ); 
		chrg_total_dist.push_back( vector< double > () ); 
	}
	
	//#pragma omp parallel for private(grid_occ)
	
	// Start the BKL algorithm
	for (int i=0; i<n_frame; i++){
		
		// Variables for each charge
		vector<int> curr_mol, curr_box; // Number of the molecule in the mini-grid, and number of this mini-grid
		vector<double> chrg_tmp_dist, chrg_tmp_jump; // Distance and number of jumps for each charge
		
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
		// Set charge variables to zero
		total_time_try = 0.0;
		total_dist_try = 0.0;
		chrg_tmp_dist.clear();
		chrg_tmp_jump.clear(); 
		
		for (unsigned int charge_i=0; charge_i<n_charges; charge_i++){
			chrg_tmp_dist.push_back(0.0);
			chrg_tmp_jump.push_back(0.0);
		
			chrg_E_electrostatic[i].push_back( vector < double > () ); 
			chrg_E_0[i].push_back( vector < double > () );
			chrg_E_1[i].push_back( vector < double > () );
			chrg_E_electrostatic[i][charge_i].push_back(Calcul_V(i, curr_mol[charge_i], charge_i, curr_mol, curr_box));
			chrg_total_time[i].push_back(0.0);
			chrg_total_dist[i].push_back(0.0);
			
			if (grid_E_type.compare(0,5,"INPUT") == 0) {
				chrg_E_0[i][charge_i].push_back(E_0[i][curr_mol[charge_i]]); 
				chrg_E_1[i][charge_i].push_back(E_1[i][curr_mol[charge_i]]);
			}
			else {
				chrg_E_0[i][charge_i].push_back(E_grid[i][curr_box[charge_i]][curr_mol[charge_i]]); 
				chrg_E_1[i][charge_i].push_back(0.0);
			}
		}

		bool previous_jump_ok = true;
		bool exit_loop = false;
		double sum_k = 0.0;
		
		// ---------------------------------------------------------------------

		unsigned int charge_try = 0;
		while (charge_try<n_try*n_charges){	
			
			// Calculates the transfer rates
			if (previous_jump_ok){
				
				// List all the possible events
				event_k.clear();
				event_charge.clear();
				event_mol_index.clear();
				event_neigh_num.clear(); 
				event_neigh_index.clear();
				
				for (unsigned int charge_i=0; charge_i<curr_mol.size(); charge_i++){
					
					int tmp_mol_index = curr_mol[charge_i];
					int tmp_mol_box = curr_box[charge_i];
											
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
						
						double k_tmp = numeric_limits<double>::max();
						double dE_tmp = 0.0;
						
						// Choose deltaE depending on the type of disorder
						if (grid_E_type.compare(0,5,"INPUT") == 0) {
							dE_tmp = dE_box[i][tmp_mol_index][jj];
						}
						
						else {
							dE_tmp = dE_grid[i][tmp_mol_box][tmp_mol_index][jj];
						}
						
						// Calculate transfer rate
						if (tmp_neigh_index != tmp_mol_index){
							k_tmp = Marcus_Levich_Jortner_rate_electro(i, tmp_mol_index, tmp_neigh_index, tmp_neigh_num, d_x[i][tmp_mol_index][jj], d_y[i][tmp_mol_index][jj], d_z[i][tmp_mol_index][jj], dE_tmp, J_H[i][tmp_mol_index][jj], J_L[i][tmp_mol_index][jj], curr_mol, curr_box, charge_i);
						}

						event_k.push_back(k_tmp);
						event_charge.push_back(charge_i);
						event_mol_index.push_back(tmp_mol_index);
						event_neigh_num.push_back(jj);
						event_neigh_index.push_back(tmp_neigh_index);
					}
				}
				
				// Sum of transfer rates of all the events
				sum_k = 0.0;
				for (unsigned int event=0; event<event_k.size(); event++){
					sum_k += event_k[event];
				}
				
				exit_loop = false;
			}
			
			else{
				exit_loop = false;
			}
			
			// Choose an event randomly
			//double random_k = Rand_0_1() * sum_k;
			double random_k = Rand_0_1_boost(rng) * sum_k;
			double partial_sum_k = 0.0;
			
			// Find this event in the list
			for (unsigned int event=0; event<event_k.size() && exit_loop == false; event++){
				partial_sum_k += event_k[event];
				
				if (partial_sum_k > random_k){
					
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
						//printf("out of system\n");
					}
					
					// Do not jump if the next molecule is occupied
					else if (grid_occ[tmp_curr_box][tmp_curr_mol] == true){
						previous_jump_ok = false;
						exit_loop = true;
						//printf("occupied\n");
					}
					
					// Do not jump if the charge change of layer
					else if (n_layer > 0 && list_layer[i][curr_mol[event_charge[event]]] != list_layer[i][tmp_curr_mol]){
						previous_jump_ok = false;
						exit_loop = true;
					}
					
					else{
						
						// Calculate the time for each charge and the total time
						//double k_event = -log(Rand_0_1())/(sum_k);
						double k_event = -log(Rand_0_1_boost(rng))/(sum_k);
						chrg_total_time[i][event_charge[event]] += k_event;
						total_time_try += k_event;
						
						// Calculate the distance traveled by the charge and the total distance
						double event_dist = (d_x[i][event_mol_index[event]][event_neigh_num[event]] * uF_x + d_y[i][event_mol_index[event]][event_neigh_num[event]] * uF_y + d_z[i][event_mol_index[event]][event_neigh_num[event]] * uF_z)*1E-8;
						chrg_total_dist[i][event_charge[event]] += event_dist/double(curr_mol.size());
						chrg_tmp_dist[event_charge[event]] += event_dist/double(curr_mol.size());
						total_dist_try += event_dist/double(curr_mol.size());
						// Note:
						// We divide here by the number of charges in the system, 
						// so there is no problem afterwards if the number of charges
						// in the system change. Indeed:
						// <x> = (1/n) * sum_i x_i = sum_i (x_i/n)
									
						// Calculate the number of jumps for each charge
						chrg_tmp_jump[event_charge[event]] += 1.0;
						
						// Set the occupancy of the grid
						grid_occ[curr_box[event_charge[event]]][event_mol_index[event]] = false;
						grid_probability[i][curr_box[event_charge[event]]][event_mol_index[event]][event_charge[event]] += event_dist;
						
						// Set the new position in the grid from temp values
						curr_box[event_charge[event]] = tmp_curr_box;
						curr_mol[event_charge[event]] = tmp_curr_mol;
						
						// Check if the charge reached the end of the grid
						if ((F_dir.compare("a") == 0 && tmp_curr_box_a >= n_box_a-1) || (F_dir.compare("b") == 0 && tmp_curr_box_b >= n_box_b-1) || (F_dir.compare("c") == 0 && tmp_curr_box_c >= n_box_c-1) || (F_dir.compare("ab") == 0 && fabs(chrg_tmp_dist[event_charge[event]]) > dist_tot) || (F_dir.compare("ac") == 0 && fabs(chrg_tmp_dist[event_charge[event]]) > dist_tot) || (F_dir.compare("bc") == 0 && fabs(chrg_tmp_dist[event_charge[event]]) > dist_tot)){							
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
							chrg_tmp_dist[event_charge[event]] = 0.0;
							
							delete [] pos;
							
							// Set distance and jumps to zero
							chrg_tmp_dist[event_charge[event]] = 0.0;
							chrg_tmp_jump[event_charge[event]] = 0.0;
							
							// Print summary for the current try in the file
							if ((charge_try/n_charges) == (double(charge_try)/double(n_charges))) {
								Print_Summary_Try(output_folder, i, charge_try, total_dist_try, total_time_try);
							}
							
							// Print info on standard output
							if ((charge_try/n_charges) % 10 == 0 && (charge_try/n_charges) == (double(charge_try)/double(n_charges))) {
								t_info = time(NULL); t_info_str = asctime(localtime(&t_info)); 
								t_info_str.erase(t_info_str.length()-1,1); 
								cout << "[INFO: " << t_info_str << "] Running Frame " << i+1 << "/" << n_frame << ", Try " << (charge_try/n_charges) << "/" << n_try << endl;
							}
							
							charge_try++;
							
						}
						
						// Set grid and charge properties if the limit is not reached
						else{
							grid_occ[tmp_curr_box][tmp_curr_mol] = true;
							chrg_E_electrostatic[i][event_charge[event]].push_back(Calcul_V(i, event_neigh_index[event], event_charge[event], curr_mol, curr_box));  
							if (grid_E_type.compare(0,5,"INPUT") == 0) {
								chrg_E_0[i][event_charge[event]].push_back(E_0[i][tmp_curr_mol]); 
								chrg_E_1[i][event_charge[event]].push_back(E_1[i][tmp_curr_mol]);
							}
							else {
								chrg_E_0[i][event_charge[event]].push_back(E_grid[i][tmp_curr_box][tmp_curr_mol]); 
								chrg_E_1[i][event_charge[event]].push_back(0.0);
							}
						}
						
						previous_jump_ok = true;
						exit_loop = true;
					}
				}
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
		chrg_tmp_dist.clear(); 
		chrg_tmp_jump.clear();
		
		event_k.clear();
		event_charge.clear();
		event_mol_index.clear();
		event_neigh_num.clear();
		event_neigh_index.clear();
	}
	
	//#pragma omp barrier
	
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
