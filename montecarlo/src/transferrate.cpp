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
#include <iostream> //Entrées-sorties standard
#include <string> //Chaines de caracteres
#include <limits> //Pour aller à la fin d'une ligne, par exemple
#include <vector>

// C libraries
#include <math.h>

// Local libraries
#include "constants.h"
#include "variables.h"
#include "mathplus.h"
#include "deltav.h"

using namespace std;

void Marcus_Levich_Jortner_CST(){
		
	S = LAMBDA_I/H_OMEGA;
	MLJ_CST1 = (2*PI/H_BAR)*(1.0/(sqrt(4*PI*LAMBDA_S*K_BOLTZ*T)));
	MLJ_CST2 = 4*LAMBDA_S*K_BOLTZ*T;
	
	int n=0;
	double fact = Facto(n);
	double tmp = exp(-S)*(pow(S,n)/fact);
	
	while (fact < numeric_limits<double>::max()) {
		MLJ_CST3.push_back(tmp);
		n++;
		fact *= n;
		tmp = exp(-S)*(pow(S,n)/(fact));
	}
}

double Marcus_Levich_Jortner_rate(double d_x_tmp, double d_y_tmp,\
				double d_z_tmp, double dE_tmp, double J_H_tmp, double J_L_tmp){
	
	double dG0 = 0.0;	
	
	// CHECK SIGN!
	if (charge.compare("e") == 0)
		dG0 = -(d_x_tmp * F_x + d_y_tmp * F_y + d_z_tmp * F_z)*1E-8;
		
	else if (charge.compare("h") == 0)
		dG0 = (d_x_tmp * F_x + d_y_tmp * F_y + d_z_tmp * F_z)*1E-8;
		
	dG0 = dG0 + dE_tmp;
	
	int n = 0;
	double k_tmp = 0.0;
	double k_inter = MLJ_CST3[n]*\
							exp(-pow(dG0 + LAMBDA_S + n*H_OMEGA,2)/(MLJ_CST2));
	
	while (k_inter > numeric_limits<double>::min()) {
		k_tmp += k_inter;
		n++;
		k_inter = MLJ_CST3[n]*\
							exp(-pow(dG0 + LAMBDA_S + n*H_OMEGA,2)/(MLJ_CST2));
	}
	
	if (charge.compare("e") == 0)
		k_tmp = MLJ_CST1 * pow(J_L_tmp,2) * k_tmp;
	
	else if (charge.compare("h") == 0)
		k_tmp = MLJ_CST1 * pow(J_H_tmp,2) * k_tmp;
		
	return k_tmp;
	
}

double Marcus_Levich_Jortner_rate_electro(int i, int mol_index_tmp,\
		int neigh_index_tmp, int neigh_num_tmp, double d_x_tmp, double d_y_tmp,\
		double d_z_tmp, double dE_tmp, double J_H_tmp, double J_L_tmp,\
		vector<int> curr_mol_tmp, vector<int> curr_box_tmp,\
		unsigned int charge_i_tmp){
	
	double dV = 0.0;
	double dG0 = 0.0;	
	
	// Calcul DeltaV
	dV = Calcul_DeltaV(i, mol_index_tmp, neigh_index_tmp, neigh_num_tmp,\
									charge_i_tmp, curr_mol_tmp, curr_box_tmp);

	// CHECK SIGN!
	if (charge.compare("e") == 0)
		dG0 = -(d_x_tmp * F_x + d_y_tmp * F_y + d_z_tmp * F_z)*1E-8;
		
	else if (charge.compare("h") == 0)
		dG0 = (d_x_tmp * F_x + d_y_tmp * F_y + d_z_tmp * F_z)*1E-8;
	
	//printf("%e %e\n", dG0, dV);	
	
	dG0 = dG0 + dE_tmp + dV;
	
	int n = 0;
	double k_tmp = 0.0;
	double k_inter = MLJ_CST3[n]*\
							exp(-pow(dG0 + LAMBDA_S + n*H_OMEGA,2)/(MLJ_CST2));
	
	while (k_inter > numeric_limits<double>::min()) {
		k_tmp += k_inter;
		n++;
		k_inter = MLJ_CST3[n]*\
							exp(-pow(dG0 + LAMBDA_S + n*H_OMEGA,2)/(MLJ_CST2));
	}
	
	if (charge.compare("e") == 0)
		k_tmp = MLJ_CST1 * pow(J_L_tmp,2) * k_tmp;
	
	else if (charge.compare("h") == 0)
		k_tmp = MLJ_CST1 * pow(J_H_tmp,2) * k_tmp;
		
	return k_tmp;
	
}

// Calculates transfer rates between molecules
void Calcul_k(bool print_results){
	
	k.clear();
	
	for (int i=0; i<n_frame; i++){
		k.push_back( vector< vector<double> > ());
		
		for (int ii=0; ii<n_mol; ii++){
			k[i].push_back( vector<double> ());
			
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
				
				double k_tmp = Marcus_Levich_Jortner_rate(d_x[i][ii][jj],d_y[i][ii][jj], d_z[i][ii][jj],\
															dE[i][ii][jj], J_H[i][ii][jj], J_L[i][ii][jj]);
				
				k[i][ii].push_back(k_tmp);
				
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
					cout << neigh_label[i][ii][jj] << " " << k[i][ii][jj] << endl;
				}
			}
		}
	}
}

// Calculates the full matrix (all neighbors of all molecules)
void Full_Matrix(){
	
	for (int i=0; i<n_frame; i++){
		for (int ii=1; ii<n_mol; ii++){
			for (int ll=ii-1; ll>=0; ll--){
				
				//Loop on all the molecules to find the missing neighbors
				for (unsigned int jj=0; jj<neigh_label[i][ll].size(); jj++){
					
					if (neigh_label[i][ll][jj]==mol_label[ii]){
						
						neigh_label[i][ii].insert(neigh_label[i][ii].begin(), mol_label[ll]);
						
						d_x[i][ii].insert(d_x[i][ii].begin(), -d_x[i][ll][jj]);
						d_y[i][ii].insert(d_y[i][ii].begin(), -d_y[i][ll][jj]);
						d_z[i][ii].insert(d_z[i][ii].begin(), -d_z[i][ll][jj]);
						
						neigh_jump_vec_a[i][ii].insert(neigh_jump_vec_a[i][ii].begin(),\
																				-neigh_jump_vec_a[i][ll][jj]);
						neigh_jump_vec_b[i][ii].insert(neigh_jump_vec_b[i][ii].begin(),\
																				-neigh_jump_vec_b[i][ll][jj]);
						neigh_jump_vec_c[i][ii].insert(neigh_jump_vec_c[i][ii].begin(),\
																				-neigh_jump_vec_c[i][ll][jj]);

						dE[i][ii].insert(dE[i][ii].begin(), -dE[i][ll][jj]);
						
						J_H[i][ii].insert(J_H[i][ii].begin(), J_H[i][ll][jj]);
						J_L[i][ii].insert(J_L[i][ii].begin(), J_L[i][ll][jj]);
							
						double k_tmp = Marcus_Levich_Jortner_rate(d_x[i][ii].front(), d_y[i][ii].front(),\
							d_z[i][ii].front(), dE[i][ii].front(), J_H[i][ii].front(), J_L[i][ii].front());
						
						k[i][ii].insert(k[i][ii].begin(), k_tmp);

						break;
					}
				}
			}
		}
	}	
}

// Calculates 1/k, and set small k to zero
void Inverse_Clear_k(bool print_results){
	
	k_inv.clear();
	
	for (int i=0; i<n_frame; i++){
		k_inv.push_back( vector< vector<double> > ());
		
		for (int ii=0; ii<n_mol; ii++){
			k_inv[i].push_back( vector<double> ());
			
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
				
				if (k[i][ii][jj] < 1E+08)
					k_inv[i][ii].push_back(numeric_limits<double>::max());
				else
					k_inv[i][ii].push_back(1.0/k[i][ii][jj]);
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
					cout << neigh_label[i][ii][jj] << " " << k_inv[i][ii][jj] << endl;
				}
			}
		}
	}
}
