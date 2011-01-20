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

#ifndef _TRANSFERRATE_H
#define _TRANSFERRATE_H 1

// Estimates some constant parts on the MLJ rate
void Marcus_Levich_Jortner_CST();

// MLJ rate calculation
double Marcus_Levich_Jortner_rate(double d_x_tmp, double d_y_tmp, double d_z_tmp, double dE_tmp,\
																			double J_H_tmp, double J_L_tmp);

// MLJ rate with electrostatic interactions calculation
double Marcus_Levich_Jortner_rate_electro(int i, int mol_index_tmp,\
						int neigh_index_tmp, int neigh_num_tmp,\
						double d_x_tmp, double d_y_tmp, double d_z_tmp,\
						double dE_tmp, double J_H_tmp, double J_L_tmp,\
						vector<int> curr_mol_tmp, vector<int> curr_box_tmp,\
						unsigned int charge_i_tmp);

// Miller Abrahams rate with electrostatic interactions calculation					
double Miller_Abrahams_rate_electro(int i, int mol_index_tmp, int neigh_index_tmp, int neigh_num_tmp,\
		double d_x_tmp, double d_y_tmp, double d_z_tmp, double dE_tmp, double J_H_tmp, double J_L_tmp,\
		vector<int> curr_mol_tmp, vector<int> curr_box_tmp, unsigned int charge_i_tmp);

// Calculates transfer rates between molecules
void Calcul_k(bool print_results);

// Calculates the full matrix (all neighbors of all molecules)
void Full_Matrix();

// Calculates 1/k, and set small k to zero
void Inverse_Clear_k(bool print_results);


#endif // transferrate.h
