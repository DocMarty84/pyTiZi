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
#include <iomanip>
#include <sstream>
#include <vector>

// C libraries
#include <stdlib.h>

// Local libraries
#include "clear.h"
#include "variables.h"

using namespace std;

// =============================================================================
// --------------------------------- Read part ---------------------------------
// =============================================================================

// Read .mc file
void Read_MC(string input_file, string input_folder, bool print_results){
	int n;
	int mol_label_tmp, neigh_label_tmp;
	double J_H_tmp, J_L_tmp;
	
	string tmp;
	stringstream file_mc;
	file_mc << input_folder << "/" << input_file << ".mc";
	
	ifstream input(file_mc.str().c_str(), ios::in);
	if (input){
		
		/* 
		 * ===============================================
		 * Physical units:
		 * ---------------
		 * 
		 * snap_delay: seconds
		 * LAMBDA_I_H, LAMBDA_I_E, LAMBDA_S, H_OMEGA: eV
		 * T: Kelvin
		 * dist_tot: centimeter
		 * F_norm: V/cm
		 * ===============================================
		 */
		
		input >> n_frame >> n_mol;
		input >> snap_delay;
		input >> LAMBDA_I_H >> LAMBDA_I_E >> LAMBDA_S >> T >> H_OMEGA >> dist_tot >> n_try;
		//input >> n_box_a >> n_box_b >> n_box_c >> n_charges;
		input >> n_box_a >> n_box_b >> n_box_c >> tmp;
		input >> F_norm;
		
		n_box = n_box_a * n_box_b * n_box_c;
		
		//Generate vectors
		for (int i=0; i<n_frame; i++){
			neigh_label.push_back( vector< vector<int> > ());
			J_H.push_back( vector< vector<double> > ());
			J_L.push_back( vector< vector<double> > ());
			
			for (int ii=0; ii<n_mol; ii++){
				neigh_label[i].push_back( vector<int> ());
				J_H[i].push_back( vector<double> ());
				J_L[i].push_back( vector<double> ());
			}
		}
		
		//Read J
		for (int i=0; i<n_frame; i++){
			input >> tmp >> tmp;
			
			for (int ii=0; ii<n_mol; ii++){
				input >> tmp >> mol_label_tmp >> n;
				mol_label.push_back(mol_label_tmp);
				
				for (int jj=0; jj<n; jj++){
					input >> neigh_label_tmp >> J_H_tmp >> J_L_tmp;
					
					neigh_label[i][ii].push_back(neigh_label_tmp);
					J_H[i][ii].push_back(J_H_tmp);
					J_L[i][ii].push_back(J_L_tmp);
				}
			}
		}
		
		input.close();
		
		// Print part
		if (print_results){
		cout << n_frame << " " << n_mol << endl;
		cout << snap_delay << endl;
		cout << LAMBDA_I << " " << LAMBDA_S << " " << T << " " << H_OMEGA << " " << dist_tot << " "\
																					<< n_try << " " << endl;
		cout << n_box_a << " " << n_box_b << " " << n_box_c << " " << n_charges << endl;
		cout << F_norm << endl;		

			for (int i=0; i<n_frame; i++){
				cout << "frame " << i << endl;
				
				for (int ii=0; ii<n_mol; ii++){
					input >> tmp >> mol_label_tmp >> n;
					mol_label.push_back(mol_label_tmp);
					cout << "molecule " << mol_label[ii] << " " << neigh_label[i][ii].size() << endl;
					
					for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
						cout << neigh_label[i][ii][jj] << " " << J_H[i][ii][jj] << " " << J_L[i][ii][jj]\
																									<< endl;
					}
				}
			}
		}
		
	}
	else{
		cerr << "Error opening " << file_mc.str().c_str() << endl;
		Clear_All();
		exit(1);
	}
}

// Read .cell file
void Read_CELL(string input_file, string input_folder, bool print_results){
	double a_tmp, b_tmp, c_tmp;
	double alpha_deg_tmp, beta_deg_tmp, gamma_deg_tmp, vol_box_tmp;
	double temp_alpha_cos_tmp, temp_beta_sin_tmp, temp_beta_cos_tmp, temp_gamma_sin_tmp, temp_gamma_cos_tmp,\
																	temp_beta_term_tmp, temp_gamma_term_tmp;

	string tmp;
	stringstream file_cell;
	file_cell << input_folder << "/" << input_file << ".cell";
	
	ifstream input(file_cell.str().c_str(), ios::in);
	if (input){

		input >> tmp >> pbc[0] >> pbc[1] >> pbc[2];
		input >> tmp >> tmp;
		for(int i=0; i<n_frame; i++){
			input >> tmp >> tmp;
			input >> a_tmp >> b_tmp >> c_tmp >> alpha_deg_tmp >> beta_deg_tmp >> gamma_deg_tmp >> vol_box_tmp;
			input >> temp_alpha_cos_tmp >> temp_beta_sin_tmp >> temp_beta_cos_tmp >> temp_gamma_sin_tmp >>\
											temp_gamma_cos_tmp >> temp_beta_term_tmp >> temp_gamma_term_tmp;
			
			a.push_back(a_tmp); 
			b.push_back(b_tmp); 
			c.push_back(c_tmp);
			alpha_deg.push_back(alpha_deg_tmp); 
			beta_deg.push_back(beta_deg_tmp); 
			gamma_deg.push_back(gamma_deg_tmp);
			vol_box.push_back(vol_box_tmp);

			temp_alpha_cos.push_back(temp_alpha_cos_tmp);
			temp_beta_sin.push_back(temp_beta_sin_tmp);
			temp_beta_cos.push_back(temp_beta_cos_tmp);
			temp_gamma_sin.push_back(temp_gamma_sin_tmp);
			temp_gamma_cos.push_back(temp_gamma_cos_tmp);
			temp_beta_term.push_back(temp_beta_term_tmp);
			temp_gamma_term.push_back(temp_gamma_term_tmp);

		}
		
		input.close();
		
		// Print part
		if (print_results){
			for(int i=0; i<n_frame; i++){
				cout << "frame " << i << endl;
				cout << a[i] << " " << b[i] << " " << c[i] << " " << alpha_deg[i] << " " << beta_deg[i] <<\
																				" " << gamma_deg[i] << endl;
				cout << temp_alpha_cos[i] << " " << temp_beta_sin[i] << " " << temp_beta_cos[i] << " "\
												<< temp_gamma_sin[i] << " " << temp_gamma_cos[i] << " "\
												<< temp_beta_term[i] << " " << temp_gamma_term[i] << endl;

			}
		}
		
	}
	else{
		cerr << "Error opening " << file_cell.str().c_str() << endl;
		Clear_All();
		exit(1);
	}
}

// Read .cm part
void Read_CM(string input_file, string input_folder, bool print_results){
	double CM_x_tmp, CM_y_tmp, CM_z_tmp;
	
	string tmp;
	stringstream file_cm;
	file_cm << input_folder << "/" << input_file << ".cm";
	
	ifstream input(file_cm.str().c_str(), ios::in);
	if (input){
		
		// Generate vectors
		for (int i=0; i<n_frame; i++) {
			CM_x.push_back( vector<double> ());
			CM_y.push_back( vector<double> ());
			CM_z.push_back( vector<double> ());
		}
		
		// Read CM
		for(int i=0; i<n_frame; i++){
			input >> tmp >> tmp;
			for(int ii=0; ii<n_mol; ii++){
				input >> tmp >> tmp >> tmp >> CM_x_tmp >> CM_y_tmp >> CM_z_tmp >> tmp;

				CM_x[i].push_back(CM_x_tmp);
				CM_y[i].push_back(CM_y_tmp);
				CM_z[i].push_back(CM_z_tmp);
			}
		}
		
		input.close();
		
		// Print part
		if(print_results){
			for(int i=0; i<n_frame; i++){
				cout << "frame " << i << endl;
				for(int ii=0; ii<n_mol; ii++){
					cout << "molecule " << mol_label[ii] << " " << CM_x[i][ii] << " " << CM_y[i][ii] <<\
																				" " << CM_z[i][ii] << endl;
				}
			}
		}
	}
	else{
		cerr << "Error opening " << file_cm.str().c_str() << endl;
		Clear_All();
		exit(1);
	}
}

// Read e_av part (average energies)
void Read_E_av(string input_file, string input_folder, bool print_results){
	
	stringstream file_e_av;
	file_e_av << input_folder << "/" << input_file << ".e_av";
	
	ifstream input(file_e_av.str().c_str(), ios::in);
	if (input){
		string tmp;
		
		input >> tmp;
		
		if (tmp.compare("INPUT_AVERAGE") == 0) {
			
			t_info = time(NULL); t_info_str = asctime(localtime(&t_info)); 
			t_info_str.erase(t_info_str.length()-1,1); 
			cout << "[INFO: " << t_info_str << "] Average energy read from file "\
																<< file_e_av.str().c_str() << " !" << endl;
			
			grid_E_type = tmp;
			double E_0_tmp, E_1_tmp;

			// Generate vectors
			for (int i=0; i<n_frame; i++) {
				E_0.push_back( vector<double> ());
				E_1.push_back( vector<double> ());
			}
			
			// Read E (same for all the frames)
			for(int ii=0; ii<n_mol; ii++){
				input >> tmp >> E_0_tmp >> E_1_tmp;
				
				for (int i=0; i<n_frame; i++) {
					E_0[i].push_back(E_0_tmp);
					E_1[i].push_back(E_1_tmp);
				}
			}
			
			input.close();
			
			// Print part
			if(print_results){
				for(int ii=0; ii<n_mol; ii++){
					cout << "molecule " << mol_label[ii] << " " << E_0[0][ii] << " " << E_1[0][ii] << endl;
				}
			}
			
		}
		
		else if (tmp.compare("GAUSSIAN") == 0) {
			
			t_info = time(NULL); t_info_str = asctime(localtime(&t_info)); 
			t_info_str.erase(t_info_str.length()-1,1); 
			cout << "[INFO: " << t_info_str << "] Gaussian DOS for the grid will be calculated!" << endl;
			
			grid_E_type = tmp;
			input >> tmp;
			grid_sigma_over_kT = atof(tmp.c_str());
			
			input.close();
		}
		
		else if (tmp.compare("CORRELATED") == 0) {
			
			t_info = time(NULL); t_info_str = asctime(localtime(&t_info)); 
			t_info_str.erase(t_info_str.length()-1,1); 
			cout << "[INFO: " << t_info_str << "] Correlated DOS for the grid will be calculated!" << endl;
			
			grid_E_type = tmp;
			input >> tmp;
			grid_sigma_over_kT = atof(tmp.c_str());
			
			input.close();
		}
		
		else if (tmp.compare("HARD_SPHERE") == 0) {
			
			t_info = time(NULL); t_info_str = asctime(localtime(&t_info)); 
			t_info_str.erase(t_info_str.length()-1,1); 
			cout << "[INFO: " << t_info_str << "] Hard Sphere DOS for the grid will be calculated!" << endl;
			
			grid_E_type = tmp;
			input >> tmp;
			grid_radius_sphere = atof(tmp.c_str());
			
			input.close();
		}
		
		else if (tmp.compare("EXPONENTIAL_SPHERE") == 0) {
			
			t_info = time(NULL); t_info_str = asctime(localtime(&t_info)); 
			t_info_str.erase(t_info_str.length()-1,1); 
			cout << "[INFO: " << t_info_str << "] Exponential Sphere DOS for the grid will be calculated!"\
																									<< endl;
			
			grid_E_type = tmp;
			input >> tmp;
			grid_e_max_sphere = atof(tmp.c_str());
			input >> tmp;
			grid_e_decrease_sphere = atof(tmp.c_str());
			
			input.close();
		}
		
		else {
			
			cout << "[WARNING] Bad keyword specified in average energy file " << file_e_av.str().c_str() <<\
																	" !\n[WARNING] E set to zero." << endl;
			
			grid_E_type = "INPUT_ZEROS";
			
			for (int i=0; i<n_frame; i++) {
				E_0.push_back( vector<double> (n_mol, 0));
				E_1.push_back( vector<double> (n_mol, 0));
			}
			
			input.close();
		}
	}

	

	else {
		cout << "[WARNING] Average energy file " << file_e_av.str().c_str() <<\
															" not found!\n[WARNING] E set to zero." << endl;
		
		for (int i=0; i<n_frame; i++) {
			E_0.push_back( vector<double> (n_mol, 0));
			E_1.push_back( vector<double> (n_mol, 0));
		}
	}
}
