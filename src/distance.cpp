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

// C libraries
#include <math.h>

// Local libraries
#include "variables.h"
#include "coordinates.h"

using namespace std;

// Calculates distance between molecules
void Calcul_Dist(bool print_results){
	double *Dist_Cart, *Dist_Frac, *vec;
	Dist_Cart = new double[3];
	Dist_Frac = new double[3];
	vec = new double[3];
	
	d_x.clear(); d_y.clear(); d_z.clear();
	neigh_jump_vec_a.clear(); neigh_jump_vec_b.clear(); neigh_jump_vec_c.clear();
	
	for (int i=0; i<n_frame; i++){
		d_x.push_back( vector< vector<double> > ());
		d_y.push_back( vector< vector<double> > ());
		d_z.push_back( vector< vector<double> > ());
		neigh_jump_vec_a.push_back( vector< vector<int> > ());
		neigh_jump_vec_b.push_back( vector< vector<int> > ());
		neigh_jump_vec_c.push_back( vector< vector<int> > ());
		
		for (int ii=0; ii<n_mol; ii++){
			d_x[i].push_back( vector<double> ());
			d_y[i].push_back( vector<double> ());
			d_z[i].push_back( vector<double> ());
			neigh_jump_vec_a[i].push_back( vector<int> ());
			neigh_jump_vec_b[i].push_back( vector<int> ());
			neigh_jump_vec_c[i].push_back( vector<int> ());
			
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
				
				//Loop on all the molecules to find the CM of the neighbor
				for (int ll=ii+1; ll<n_mol; ll++){
					if (mol_label[ll]==neigh_label[i][ii][jj]){
						
						//Distance calculation
						Dist_Cart[0] = CM_x[i][ll] - CM_x[i][ii];
						Dist_Cart[1] = CM_y[i][ll] - CM_y[i][ii];
						Dist_Cart[2] = CM_z[i][ll] - CM_z[i][ii];
						
						Cartesian_To_Fractional(Dist_Cart, Dist_Frac, i);
					
						for (int k=0; k<3; k++){
							vec[k] = 0.0;
							
							if (fabs(Dist_Frac[k]) > 0.5 && pbc[k]){
								//cout << "check mol " << mol_label[ii] << " " << mol_label[ll] << endl;
								if (Dist_Frac[k] < 0.0)
									vec[k] = 1.0;
								else
									vec[k] = -1.0;
							}
							
							Dist_Frac[k] = Dist_Frac[k] + vec[k];
							
						}
						
						Fractional_To_Cartesian(Dist_Frac, Dist_Cart, i);
						
						d_x[i][ii].push_back(Dist_Cart[0]);
						d_y[i][ii].push_back(Dist_Cart[1]);
						d_z[i][ii].push_back(Dist_Cart[2]);
						neigh_jump_vec_a[i][ii].push_back(int(vec[0]));
						neigh_jump_vec_b[i][ii].push_back(int(vec[1]));
						neigh_jump_vec_c[i][ii].push_back(int(vec[2]));
						
						break;
					}
				}
			}
		}
	}
	
	delete [] Dist_Cart;
	delete [] Dist_Frac;
	delete [] vec;
	
	// Print part
	if (print_results){
		
		for (int i=0; i<n_frame; i++){
			cout << "frame " << i << endl;
			for (int ii=0; ii<n_mol; ii++){
				cout << "molecule " << mol_label[ii] << endl;
				for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
					double d = sqrt(pow(d_x[i][ii][jj],2) + pow(d_y[i][ii][jj],2) + pow(d_z[i][ii][jj],2));
					cout << neigh_label[i][ii][jj] << " " << d_x[i][ii][jj] << " " << d_y[i][ii][jj] << " " << d_z[i][ii][jj] << " " << d << " " << neigh_jump_vec_a[i][ii][jj] << " " << neigh_jump_vec_b[i][ii][jj] << " " << neigh_jump_vec_c[i][ii][jj] << endl;
				}
			}
		}
	}
}
