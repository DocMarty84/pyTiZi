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
#include <iostream>
#include <numeric> 
#include <vector>

// C libraries
#include <math.h>

#include <omp.h>

// Local libraries
#include "coordinates.h"
#include "constants.h"
#include "variables.h"

using namespace std;

// Calculates deltaV (electrostatic interactions) between molecules
double Calcul_DeltaV(int i, int mol_index_tmp, int neigh_index_tmp, int neigh_num_tmp,\
							unsigned int charge_i_tmp, vector<int> curr_mol_tmp, vector<int> curr_box_tmp){
	
	int max_thr = min(int(n_charges), int(omp_get_max_threads()));
	
	const double CST1 = 1.0/(4.0*PI*EPSILON_0*EPSILON_R);
	const double CUTOFF_2 = pow(CUTOFF_ELECTRO,2); //Cutoff electro square
	
	double V_electro;
	vector<double> V_electro_thr(max_thr, 0.0);
	
	vector<double> CM_1_Cart(3, 0.0), CM_1_Frac(3, 0.0), CM_1_Neigh_Cart(3, 0.0), CM_1_Neigh_Frac(3, 0.0);
	
	CM_1_Cart[0] = CM_x[i][mol_index_tmp];
	CM_1_Cart[1] = CM_y[i][mol_index_tmp];
	CM_1_Cart[2] = CM_z[i][mol_index_tmp];
	Cartesian_To_Fractional(CM_1_Cart, CM_1_Frac, i);
	
	CM_1_Frac[0] += double(box_a[curr_box_tmp[charge_i_tmp]]);
	CM_1_Frac[1] += double(box_b[curr_box_tmp[charge_i_tmp]]);
	CM_1_Frac[2] += double(box_c[curr_box_tmp[charge_i_tmp]]);
	
	CM_1_Neigh_Cart[0] = CM_x[i][neigh_index_tmp];
	CM_1_Neigh_Cart[1] = CM_y[i][neigh_index_tmp];
	CM_1_Neigh_Cart[2] = CM_z[i][neigh_index_tmp];
	Cartesian_To_Fractional(CM_1_Neigh_Cart, CM_1_Neigh_Frac, i);
	
	CM_1_Neigh_Frac[0] += double(box_a[curr_box_tmp[charge_i_tmp]]) + neigh_jump_vec_a[i][mol_index_tmp][neigh_num_tmp];
	CM_1_Neigh_Frac[1] += double(box_b[curr_box_tmp[charge_i_tmp]]) + neigh_jump_vec_b[i][mol_index_tmp][neigh_num_tmp];
	CM_1_Neigh_Frac[2] += double(box_c[curr_box_tmp[charge_i_tmp]]) + neigh_jump_vec_c[i][mol_index_tmp][neigh_num_tmp];
	
	#pragma omp parallel for num_threads(max_thr)
	
	for (int charge_i = 0; charge_i < int(curr_mol_tmp.size()); charge_i++){
		
		// Variable definition
		double dist, dist_neigh, dist_2, dist_neigh_2;
		double V_mol = 0.0, V_neigh = 0.0;
		
		vector<double> CM_2_Cart(3, 0.0), CM_2_Frac(3, 0.0);
		vector<double> Dist_Cart(3, 0.0), Dist_Frac(3, 0.0);
		
		if (charge_i != int(charge_i_tmp)){
			CM_2_Cart[0] = CM_x[i][curr_mol_tmp[charge_i]];
			CM_2_Cart[1] = CM_y[i][curr_mol_tmp[charge_i]];
			CM_2_Cart[2] = CM_z[i][curr_mol_tmp[charge_i]];
			Cartesian_To_Fractional(CM_2_Cart, CM_2_Frac, i);
			
			CM_2_Frac[0] += double(box_a[curr_box_tmp[charge_i]]);
			CM_2_Frac[1] += double(box_b[curr_box_tmp[charge_i]]);
			CM_2_Frac[2] += double(box_c[curr_box_tmp[charge_i]]);
			
			Dist_Frac[0] = CM_2_Frac[0] - CM_1_Frac[0];
			Dist_Frac[1] = CM_2_Frac[1] - CM_1_Frac[1];
			Dist_Frac[2] = CM_2_Frac[2] - CM_1_Frac[2];
			Fractional_To_Cartesian(Dist_Frac, Dist_Cart, i);
			
			dist_2 = pow(Dist_Cart[0],2) + pow(Dist_Cart[1],2) + pow(Dist_Cart[2],2);
			dist = sqrt(dist_2);
			
			if (dist_2 > CUTOFF_2){
				
				V_mol = CST1/dist;
				
				Dist_Frac[0] = CM_2_Frac[0] - CM_1_Neigh_Frac[0];
				Dist_Frac[1] = CM_2_Frac[1] - CM_1_Neigh_Frac[1];
				Dist_Frac[2] = CM_2_Frac[2] - CM_1_Neigh_Frac[2];
				Fractional_To_Cartesian(Dist_Frac, Dist_Cart, i);
			
				dist_neigh_2 = pow(Dist_Cart[0],2) + pow(Dist_Cart[1],2) + pow(Dist_Cart[2],2);
				dist_neigh = sqrt(dist_neigh_2);
				
				V_neigh = CST1/dist_neigh; 
				
				V_electro_thr[int(omp_get_thread_num())] += V_neigh - V_mol;
			}
		}
		
		CM_2_Cart.clear();
		CM_2_Frac.clear();
		Dist_Cart.clear();
		Dist_Frac.clear();
	}
	
	#pragma omp barrier
	
	V_electro = accumulate(V_electro_thr.begin(), V_electro_thr.end(), 0.0);
	
	CM_1_Cart.clear();
	CM_1_Frac.clear();
	CM_1_Neigh_Cart.clear();
	CM_1_Neigh_Frac.clear();
	V_electro_thr.clear();
	
	//cout << V_electro << endl;
	
	return V_electro;
}

// Calculates V (electrostatic interactions) for a charge
double Calcul_V(int i, int mol_index_tmp, unsigned int charge_i_tmp, vector<int> curr_mol_tmp,\
																					vector<int> curr_box_tmp){
	
	int max_thr = min(int(n_charges), int(omp_get_max_threads()));
	
	const double CST1 = 1.0/(4.0*PI*EPSILON_0*EPSILON_R);
	const double CUTOFF_2 = pow(CUTOFF_ELECTRO,2); //Cutoff electro square
	
	double V_mol;
	vector<double> V_mol_thr(max_thr, 0.0);
	
	vector<double> CM_1_Cart(3, 0.0), CM_1_Frac(3, 0.0), CM_2_Cart(3, 0.0), CM_2_Frac(3, 0.0), Dist_Cart(3, 0.0), Dist_Frac(3, 0.0);

	CM_1_Cart[0] = CM_x[i][mol_index_tmp];
	CM_1_Cart[1] = CM_y[i][mol_index_tmp];
	CM_1_Cart[2] = CM_z[i][mol_index_tmp];
	Cartesian_To_Fractional(CM_1_Cart, CM_1_Frac, i);
	
	CM_1_Frac[0] += double(box_a[curr_box_tmp[charge_i_tmp]]);
	CM_1_Frac[1] += double(box_b[curr_box_tmp[charge_i_tmp]]);
	CM_1_Frac[2] += double(box_c[curr_box_tmp[charge_i_tmp]]);
	
	#pragma omp parallel for num_threads(max_thr)
	
	for (unsigned int charge_i = 0; charge_i < curr_mol_tmp.size(); charge_i++){
		
		// Variable definition
		double dist, dist_2;
		
		vector<double> CM_2_Cart(3, 0.0), CM_2_Frac(3, 0.0);
		vector<double> Dist_Cart(3, 0.0), Dist_Frac(3, 0.0);
		
		if (charge_i != charge_i_tmp){
			CM_2_Cart[0] = CM_x[i][curr_mol_tmp[charge_i]];
			CM_2_Cart[1] = CM_y[i][curr_mol_tmp[charge_i]];
			CM_2_Cart[2] = CM_z[i][curr_mol_tmp[charge_i]];
			Cartesian_To_Fractional(CM_2_Cart, CM_2_Frac, i);
			
			CM_2_Frac[0] += double(box_a[curr_box_tmp[charge_i]]);
			CM_2_Frac[1] += double(box_b[curr_box_tmp[charge_i]]);
			CM_2_Frac[2] += double(box_c[curr_box_tmp[charge_i]]);
			
			Dist_Frac[0] = CM_2_Frac[0] - CM_1_Frac[0];
			Dist_Frac[1] = CM_2_Frac[1] - CM_1_Frac[1];
			Dist_Frac[2] = CM_2_Frac[2] - CM_1_Frac[2];
			Fractional_To_Cartesian(Dist_Frac, Dist_Cart, i);
			
			dist_2 = pow(Dist_Cart[0],2) + pow(Dist_Cart[1],2) + pow(Dist_Cart[2],2);
			dist = sqrt(dist_2);
			
			if (dist_2 > CUTOFF_2){
				
				V_mol_thr[int(omp_get_thread_num())] += CST1/dist;
			
			}
		}
		
		CM_2_Cart.clear();
		CM_2_Frac.clear();
		Dist_Cart.clear();
		Dist_Frac.clear();
	}
	
	#pragma omp barrier
	
	V_mol = accumulate(V_mol_thr.begin(), V_mol_thr.end(), 0.0);
	
	CM_1_Cart.clear();
	CM_1_Frac.clear();
	
	return V_mol;
}
