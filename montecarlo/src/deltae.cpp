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
#include <vector>
#include <boost/random.hpp>
#include <ctime>

// Local libraries
#include "constants.h"
#include "variables.h"

using namespace std;

// Calculates deltaE between molecules
void Calcul_DeltaE(bool print_results){
	
	dE.clear();
	
	for (int i=0; i<n_frame; i++){
		dE.push_back( vector< vector<double> > ());
		
		for (int ii=0; ii<n_mol; ii++){
			dE[i].push_back( vector<double> ());
			
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){

				//Loop on all the molecules to find the CM of the neighbor
				for (int ll=ii+1; ll<n_mol; ll++){
					if (mol_label[ll]==neigh_label[i][ii][jj]){
						
						// CHECK!!!
						dE[i][ii].push_back(E_1[0][ii]+E_0[0][ll] -(E_0[0][ii]+E_1[0][ll]));
						//cout << "[WARNING] The energies are supposed to be in kcal/mol, not eV!!!" << endl;
						//dE[i][ii].push_back((E_1[0][ll]+E_0[0][ii] -(E_0[0][ll]+E_1[0][ii]))/23.06056);
						
						break;
					}
				}
			}
		}
	}
	
	// Print part
	if (print_results){
		
		for (int i=0; i<n_frame; i++){
			cout << "frame " << i << endl;
			for (int ii=0; ii<n_mol; ii++){
				cout << "molecule " << mol_label[ii] << endl;
				for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
					cout << neigh_label[i][ii][jj] << " " << dE[i][ii][jj] << endl;
				}
			}
		}
	}
}

// Generate a random Energy mapping
void Generate_E_GDM(double mean, double sigma, bool print_results){

	using namespace boost;

	// Create a Mersenne twister random number generator
	// that is seeded once with #seconds since 1970
	static mt19937 rng(static_cast<unsigned> (time(0)));
 
	// Select Gaussian probability distribution
	normal_distribution<double> norm_dist(mean, sigma);
 
	// Bind random number generator to distribution, forming a function
	variate_generator<mt19937&, normal_distribution<double> > norm_dist_sampler(rng, norm_dist);

	for (int i=0; i<n_frame; i++){
		E_random.push_back( vector< vector<double> > ());
		
		for (int x=0; x<n_box; x++){
			E_random[i].push_back( vector<double> ());
		
			for (int ii=0; ii<n_mol; ii++){
				E_random[i][x].push_back( norm_dist_sampler() );
				
			}
		}
	}
	
	// Print part
	if (print_results){
		
		for (int i=0; i<n_frame; i++){
			cout << "frame " << i << endl;
			for (int x=0; x<n_box; x++){
				cout << "box " << x << endl;
				for (int ii=0; ii<n_mol; ii++){
					cout << "molecule " << mol_label[ii] << " " << E_random[i][x][ii] << endl;

				}
			}
		}
	}
}

// Calculates deltaE between molecules for random energy
void Calcul_DeltaE_GDM(bool print_results){

	Generate_E_GDM(0.0, grid_sigma_over_kT*K_BOLTZ*T, false);
	
	int neigh_box = 0, neigh_mol = 0;
	dE_random.clear();
	
	for (int i=0; i<n_frame; i++){
		dE_random.push_back( vector< vector< vector<double> > >());
		
		for (int x=0; x<n_box; x++){
			dE_random[i].push_back( vector< vector<double> >());
		
			for (int ii=0; ii<n_mol; ii++){
				dE_random[i][x].push_back( vector<double> ());
				
				for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){

					//Loop on all the molecules to find the CM of the neighbor
					for (int ll=ii+1; ll<n_mol; ll++){
						if (mol_label[ll]==neigh_label[i][ii][jj]){
							
							neigh_mol = mol_label[ll];
							// Find the box of the neighbor
							if (neigh_jump_vec_a[i][ii][jj] == 0 && neigh_jump_vec_b[i][ii][jj] == 0 && neigh_jump_vec_c[i][ii][jj] == 0){
								neigh_box = x;
							}
							
							else {
								for (unsigned int xx=0; xx<box_neigh_label[x].size(); xx++){
									if (neigh_jump_vec_a[i][ii][jj] == box_neigh_a[x][xx] && neigh_jump_vec_b[i][ii][jj] == box_neigh_b[x][xx] && neigh_jump_vec_c[i][ii][jj] == box_neigh_c[x][xx]){
										neigh_box = box_neigh_label[x][xx];
										break;
									}
								}
							}
							
							// CHECK!!!
							dE_random[i][x][ii].push_back(E_random[i][neigh_box][neigh_mol] - E_random[i][x][ii]);
							//cout << "[WARNING] The energies are supposed to be in kcal/mol, not eV!!!" << endl;
							//dE[i][ii].push_back((E_1[0][ll]+E_0[0][ii] -(E_0[0][ll]+E_1[0][ii]))/23.06056);
							
							break;
						}
					}
				}
			}
		}
	}
	
	// Print part
	if (print_results){
		
		for (int i=0; i<n_frame; i++){
			cout << "frame " << i << endl;
			for (int x=0; x<n_box; x++){
				cout << "box " << x << endl;
				for (int ii=0; ii<n_mol; ii++){
					cout << "molecule " << mol_label[ii] << endl;
					for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
						cout << neigh_label[i][ii][jj] << " " << dE_random[i][x][ii][jj] << endl;
					}
				}
			}
		}
	}
}
