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
#include <algorithm>
#include <ctime>
#include <iostream>
#include <vector>
#include <limits>
#include <boost/random.hpp>
#include <ctime>

// C libraries
#include <math.h>

// Local libraries
#include "constants.h"
#include "variables.h"
#include "coordinates.h"

using namespace std;

// Generate a Gaussian DOS
void Generate_Gaussian_DOS(double mean, double sigma, bool print_results){
	
	t_info = time(NULL); t_info_str = asctime(localtime(&t_info)); 
	t_info_str.erase(t_info_str.length()-1,1); 
	cout << "[INFO: " << t_info_str << "] Using a Gaussian DOS, with a value of sigma/kT = "\
																				<< grid_sigma_over_kT << endl;

	using namespace boost;
	
	// ---------------------------------------------------------------------

	// Create a Mersenne twister random number generator
	// that is seeded once with #seconds since 1970
	static mt19937 rng(static_cast<unsigned> (time(0)));
 
	// Select Gaussian probability distribution
	normal_distribution<double> norm_dist(mean, sigma);
 
	// Bind random number generator to distribution, forming a function
	variate_generator<mt19937&, normal_distribution<double> > norm_dist_sampler(rng, norm_dist);
	
	// ---------------------------------------------------------------------

	for (int i=0; i<n_frame; i++){
		E_grid.push_back( vector< vector<double> > ());
		
		for (int x=0; x<n_box; x++){
			E_grid[i].push_back( vector<double> ());
		
			for (int ii=0; ii<n_mol; ii++){
				E_grid[i][x].push_back( norm_dist_sampler() );
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
					cout << "molecule " << mol_label[ii] << " " << E_grid[i][x][ii] << endl;

				}
			}
		}
	}
}

// Generate a Correlated DOS
void Generate_Correlated_DOS(double sigma, bool print_results){
	
	t_info = time(NULL); t_info_str = asctime(localtime(&t_info)); 
	t_info_str.erase(t_info_str.length()-1,1); 
	cout << "[INFO: " << t_info_str << "] Using a Correlated DOS, with a value of sigma_p/kT = " << grid_sigma_over_kT << endl;

	using namespace boost;

	// ---------------------------------------------------------------------
	
	// Create a Mersenne twister random number generator
	// that is seeded once with #seconds since 1970
	static mt19937 rng(static_cast<unsigned> (time(0)));
 
	// Select uniform probability distribution
	uniform_real<> uni_11_dist(-1,1);
 
	// Bind random number generator to distribution, forming a function
	variate_generator<mt19937&, uniform_real<> > uni_11_sampler(rng, uni_11_dist);
	
	// ---------------------------------------------------------------------

	// Calculate average distance
	double d_av = 0.0;
	double d_num = 0.0;
		
	for (int i=0; i<n_frame; i++){
		for (int ii=0; ii<n_mol; ii++){
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
				d_av += sqrt(pow(d_x[i][ii][jj],2) + pow(d_y[i][ii][jj],2) + pow(d_z[i][ii][jj],2));
				d_num += 1.0;
			}
		}
	}
	d_av = d_av/d_num;

	// Generate random dipoles
	
	vector< vector< vector<double> > > p_grid_x;
	vector< vector< vector<double> > > p_grid_y;
	vector< vector< vector<double> > > p_grid_z;
	double p_norm_tmp, p_norm;
	
	// See Phys. Rev. Lett., 1998, 81, 4472 for calculating dipole norm from sigma
	p_norm = sigma*EPSILON_0*EPSILON_R*pow(d_av, 2)/(2.35);
	
	for (int i=0; i<n_frame; i++){
		p_grid_x.push_back( vector< vector<double> > () );
		p_grid_y.push_back( vector< vector<double> > () );
		p_grid_z.push_back( vector< vector<double> > () );
		
		for (int x=0; x<n_box; x++){
			p_grid_x[i].push_back( vector<double> () );
			p_grid_y[i].push_back( vector<double> () );
			p_grid_z[i].push_back( vector<double> () );
		
			for (int ii=0; ii<n_mol; ii++){
				p_grid_x[i][x].push_back( uni_11_sampler() );
				p_grid_y[i][x].push_back( uni_11_sampler() );
				p_grid_z[i][x].push_back( uni_11_sampler() );
				
				p_norm_tmp = sqrt(pow(p_grid_x[i][x].back(), 2) + pow(p_grid_y[i][x].back(), 2) +\
																			pow(p_grid_z[i][x].back(), 2));
																			
				p_grid_x[i][x].back() = (p_grid_x[i][x].back()/p_norm_tmp) * p_norm;
				p_grid_y[i][x].back() = (p_grid_y[i][x].back()/p_norm_tmp) * p_norm;
				p_grid_z[i][x].back() = (p_grid_z[i][x].back()/p_norm_tmp) * p_norm;
				
			}
		}
	}
		
	// Fill the energy grid

	for (int i=0; i<n_frame; i++){
		E_grid.push_back( vector< vector<double> > ());
		
		for (int x=0; x<n_box; x++){
			E_grid[i].push_back( vector<double> ());
		
			for (int ii=0; ii<n_mol; ii++){
				E_grid[i][x].push_back(0.0);
			}
		}
	}
	
	for (int i=0; i<n_frame; i++){
		
		#pragma omp parallel for

		for (int x=0; x<n_box; x++){
			
			vector<double> CM_1_Cart(3, 0.0), CM_1_Frac(3, 0.0), CM_2_Cart(3, 0.0), CM_2_Frac(3, 0.0);
			vector<double> Dist_Cart(3, 0.0), Dist_Frac(3, 0.0);
			double dist, dist_2;
		
			for (int ii=0; ii<n_mol; ii++){
				
				CM_1_Cart[0] = CM_x[i][ii];
				CM_1_Cart[1] = CM_y[i][ii];
				CM_1_Cart[2] = CM_z[i][ii];
				Cartesian_To_Fractional(CM_1_Cart, CM_1_Frac, i);
				
				CM_1_Frac[0] += double(box_a[x]);
				CM_1_Frac[1] += double(box_b[x]);
				CM_1_Frac[2] += double(box_c[x]);
				
				for (int y=0; y<n_box; y++){
					
					for (int ll=0; ll<n_mol; ll++){
						
						if (y != x && ll != ii) {
							
							CM_2_Cart[0] = CM_x[i][ll];
							CM_2_Cart[1] = CM_y[i][ll];
							CM_2_Cart[2] = CM_z[i][ll];
							Cartesian_To_Fractional(CM_2_Cart, CM_2_Frac, i);
							
							CM_2_Frac[0] += double(box_a[y]);
							CM_2_Frac[1] += double(box_b[y]);
							CM_2_Frac[2] += double(box_c[y]);
							
							Dist_Frac[0] = CM_2_Frac[0] - CM_1_Frac[0];
							Dist_Frac[1] = CM_2_Frac[1] - CM_1_Frac[1];
							Dist_Frac[2] = CM_2_Frac[2] - CM_1_Frac[2];
							Fractional_To_Cartesian(Dist_Frac, Dist_Cart, i);
							
							dist_2 = pow(Dist_Cart[0],2) + pow(Dist_Cart[1],2) + pow(Dist_Cart[2],2);
							dist = sqrt(dist_2);
							
							E_grid[i][x][ii] += (Dist_Cart[0]*p_grid_x[i][y][ll] +\
										Dist_Cart[1]*p_grid_y[i][y][ll] +\
										Dist_Cart[2]*p_grid_z[i][y][ll])/(EPSILON_0*EPSILON_R*pow(dist,3));
						}
					}
				}
				
				E_grid[i][x][ii] = - E_grid[i][x][ii];
			}

			CM_1_Cart.clear(); CM_1_Frac.clear();
			CM_2_Cart.clear(); CM_2_Frac.clear();
			Dist_Cart.clear(); Dist_Frac.clear();
		}

		#pragma omp barrier
	}
		
	p_grid_x.clear(); p_grid_y.clear(); p_grid_z.clear();
	
	// Print part
	if (print_results){
		
		for (int i=0; i<n_frame; i++){
			cout << "frame " << i << endl;
			for (int x=0; x<n_box; x++){
				cout << "box " << x << endl;
				for (int ii=0; ii<n_mol; ii++){
					cout << "molecule " << mol_label[ii] << " " << E_grid[i][x][ii] << endl;

				}
			}
		}
	}
}

// Generate a hard sphere DOS
void Generate_Hard_Sphere_DOS(double radius, bool print_results){
	
	t_info = time(NULL); t_info_str = asctime(localtime(&t_info)); 
	t_info_str.erase(t_info_str.length()-1,1); 
	cout << "[INFO: " << t_info_str << "] Using a Hard Sphere DOS, with a radius = "\
																			<< radius << "Angstrom" << endl;

	
	// Find center of the grid
	
	vector< vector<double> > list_x, list_y, list_z;
	vector<double> median_x, median_y, median_z;
	
	for (int i=0; i<n_frame; i++){
		list_x.push_back(vector<double> ());
		list_y.push_back(vector<double> ());
		list_z.push_back(vector<double> ());
		
		for (int x=0; x<n_box; x++){
			
			for (int ii=0; ii<n_mol; ii++){
				
				list_x[i].push_back(grid_x[i][x][ii]);
				list_y[i].push_back(grid_y[i][x][ii]);
				list_z[i].push_back(grid_z[i][x][ii]);
				
			}
		}

		sort(list_x[i].begin(), list_x[i].end());
		sort(list_y[i].begin(), list_y[i].end());
		sort(list_z[i].begin(), list_z[i].end());
		
		median_x.push_back((list_x[i].back() - list_x[i].front())/2.0);
		median_y.push_back((list_y[i].back() - list_y[i].front())/2.0);
		median_z.push_back((list_z[i].back() - list_z[i].front())/2.0);
	}
	

	for (int i=0; i<n_frame; i++){
		E_grid.push_back( vector< vector<double> > ());
		
		for (int x=0; x<n_box; x++){
			E_grid[i].push_back( vector<double> ());
		
			for (int ii=0; ii<n_mol; ii++){
				double dist, dist_2;
				dist_2 = pow(grid_x[i][x][ii]-median_x[i],2) + pow(grid_y[i][x][ii]-median_y[i],2) + pow(grid_z[i][x][ii]-median_z[i],2);
				
				if (dist_2 > pow(radius,2)) {
					E_grid[i][x].push_back( 0.0 );
				}
				else {
					E_grid[i][x].push_back( numeric_limits<double>::max() );
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
					cout << "molecule " << mol_label[ii] << " " << E_grid[i][x][ii] << endl;

				}
			}
		}
	}
}

// Set the variable dE to zero
void DeltaE_ZERO(bool print_results){
	
	dE_box.clear();
	
	for (int i=0; i<n_frame; i++){
		dE_box.push_back( vector< vector<double> > ());
		
		for (int ii=0; ii<n_mol; ii++){
			dE_box[i].push_back( vector<double> ());
			
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
				dE_box[i][ii].push_back(0.0);
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
					cout << neigh_label[i][ii][jj] << " " << dE_box[i][ii][jj] << endl;
				}
			}
		}
	}
}

// Calculates deltaE between molecules for a distribution in the box
void Calcul_DeltaE(bool print_results){
	
	dE_box.clear();
	
	for (int i=0; i<n_frame; i++){
		dE_box.push_back( vector< vector<double> > ());
		
		for (int ii=0; ii<n_mol; ii++){
			dE_box[i].push_back( vector<double> ());
			
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){

				//Loop on all the molecules to find the CM of the neighbor
				for (int ll=ii+1; ll<n_mol; ll++){
					if (mol_label[ll]==neigh_label[i][ii][jj]){
						
						// CHECK!!!
						dE_box[i][ii].push_back(E_1[0][ii]+E_0[0][ll] -(E_0[0][ii]+E_1[0][ll]));
						//cout << "[WARNING] The energies are supposed to be in kcal/mol, not eV!!!" << endl;
						//dE_box[i][ii].push_back((E_1[0][ll]+E_0[0][ii] -(E_0[0][ll]+E_1[0][ii]))/23.06056);
						
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
					cout << neigh_label[i][ii][jj] << " " << dE_box[i][ii][jj] << endl;
				}
			}
		}
	}
}

// Calculates deltaE between molecules for a distribution in the grid
void Calcul_DeltaE_GRID(bool print_results){
	
	vector< vector< vector< vector<int> > > > dE_grid_box_list;
	int neigh_box = 0, neigh_mol = 0;
	
	dE_grid.clear();
	
	for (int i=0; i<n_frame; i++){
		dE_grid.push_back( vector< vector< vector<double> > >());
		dE_grid_box_list.push_back( vector< vector< vector<int> > >());
		
		for (int x=0; x<n_box; x++){
			dE_grid[i].push_back( vector< vector<double> >());
			dE_grid_box_list[i].push_back( vector< vector<int> >());
		
			for (int ii=0; ii<n_mol; ii++){
				dE_grid[i][x].push_back( vector<double> ());
				dE_grid_box_list[i][x].push_back( vector<int> ());

				for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){

					//Loop on all the molecules to find the CM of the neighbor
					for (int ll=0; ll<n_mol; ll++){
						if (mol_label[ll]==neigh_label[i][ii][jj]){

							neigh_mol = mol_label[ll];
							bool out_of_system = true;
							// Find the box of the neighbor
							if (neigh_jump_vec_a[i][ii][jj] == 0 && neigh_jump_vec_b[i][ii][jj] == 0 &&\
																			neigh_jump_vec_c[i][ii][jj] == 0){
								neigh_box = x;
								out_of_system = false;
							}
							
							else {
								for (unsigned int xx=0; xx<box_neigh_label[x].size(); xx++){
									if (neigh_jump_vec_a[i][ii][jj] == box_neigh_a[x][xx] &&\
														neigh_jump_vec_b[i][ii][jj] == box_neigh_b[x][xx] &&\
														neigh_jump_vec_c[i][ii][jj] == box_neigh_c[x][xx]){
										neigh_box = box_neigh_label[x][xx];
										out_of_system = false;
										break;
									}
								}
							}

							
							if (out_of_system){
								dE_grid[i][x][ii].push_back(numeric_limits<double>::max());
								dE_grid_box_list[i][x][ii].push_back(numeric_limits<int>::max());
							}
							else{
								// CHECK!!!
								dE_grid[i][x][ii].push_back(E_grid[i][neigh_box][neigh_mol] -\
																							E_grid[i][x][ii]);
								dE_grid_box_list[i][x][ii].push_back(neigh_box);
							}
							//cout << "[WARNING] The energies are supposed to be in kcal/mol, not eV!!!" << endl;
							//dE_box[i][ii].push_back((E_1[0][ll]+E_0[0][ii] -(E_0[0][ll]+E_1[0][ii]))/23.06056);
							
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
						if (dE_grid_box_list[i][x][ii][jj] != numeric_limits<int>::max()) {
							cout << dE_grid_box_list[i][x][ii][jj] << " " << neigh_label[i][ii][jj] << " "\
								<< E_grid[i][x][ii] << " "\
								<< E_grid[i][dE_grid_box_list[i][x][ii][jj]][neigh_label[i][ii][jj]] << " "\
								<< dE_grid[i][x][ii][jj] << endl;
						}
					}
				}
			}
		}
	}
	
	dE_grid_box_list.clear();
}

