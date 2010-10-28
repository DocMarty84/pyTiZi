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
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <omp.h>

using namespace std;

// Various constants
const double K_BOLTZ = 8.617343e-5; //Boltzmann constant in eV
const double H_BAR = 6.58211899e-16; // h/2PI in eV
const double EPSILON_0 = (0.5)*(137.035999084)*(1.0/4.13566733e-15)*(1.0/299792458e10) ; //Epsilon_0 (Vacuum permittivity in e/(V.Ang))
const double PI = 3.14159265358979323846;
const double EPSILON_R = 1.0; //Epsilon_r (Relative permittivity in A.s/(V.cm))
const double CUTOFF_ELECTRO = 150; //Cutoff for electrostatic interactions in Angstrom

// Variables read in the input file
int n_frame, n_mol; // Number of frame and molecule
double snap_delay;
double LAMBDA_I, LAMBDA_I_H, LAMBDA_I_E, LAMBDA_S, T, H_OMEGA, dist_tot;
unsigned int n_try, n_charges;
int n_mini_grid_a, n_mini_grid_b, n_mini_grid_c;
double F_norm; 
string F_dir;
string charge;

// Variables for the unit cell parameters
bool pbc[3]; // Periodic boundary conditions
vector<double> a, b, c, alpha_deg, beta_deg, gamma_deg, vol_box; // Cell parameters
vector<double> temp_alpha_cos, temp_beta_sin, temp_beta_cos, temp_gamma_sin, temp_gamma_cos, temp_beta_term, temp_gamma_term; // Parameters for fractional coordinates

// Variables for each molecule
vector<int> mol_label;
vector< vector<double> > CM_x, CM_y, CM_z; // Center of masses
vector< vector<double> > E_0, E_1;

// Variables for the grid
vector<int> box_a, box_b, box_c; // Position of the mini-grids
vector< vector<bool> > grid_occ;
vector< vector<int> > box_neigh_a, box_neigh_b, box_neigh_c, box_neigh_label; // Neighbor of each mini-grid

// Variables for neighbors
vector< vector< vector<int> > > neigh_label;
vector< vector< vector<double> > > d_x, d_y, d_z;
vector< vector< vector<double> > > dE;
vector< vector< vector<double> > > J_H, J_L;
vector< vector< vector<int> > > neigh_jump_vec_a, neigh_jump_vec_b, neigh_jump_vec_c; // Vector for the change in mini-grid
vector< vector< vector<double> > > k, k_inv;

// Variables for the electric field direction
double theta_deg, phi_deg, theta_rad, phi_rad;
double F_x, F_y, F_z;

// Constants for the MLJ theory
double S, MLJ_CST1, MLJ_CST2; 
vector <double> MLJ_CST3;

// Variables for MT
bool MT = 1;
vector<double> event_k_1, event_k_2; // Transfer rate
vector<int> event_charge_1, event_mol_index_1, event_neigh_num_1, event_neigh_index_1, event_charge_2, event_mol_index_2, event_neigh_num_2, event_neigh_index_2; // Initial molecule and neighbor


// =============================================================================
// ------------------------ Coordinates transformations ------------------------
// =============================================================================

void Cartesian_To_Fractional(double* Dist_Cart, double* Dist_Frac, int i){
	double x = Dist_Cart[0];
	double y = Dist_Cart[1];
	double z = Dist_Cart[2];
	
	z = (z/temp_gamma_term[i]) / c[i];
	y = ((y-z*c[i]*temp_beta_term[i])/temp_gamma_sin[i]) / b[i];
	x = (x-y*b[i]*temp_gamma_cos[i]-z*c[i]*temp_beta_cos[i]) / a[i];
	Dist_Frac[0] = x; Dist_Frac[1] = y; Dist_Frac[2] = z;	
}

void Fractional_To_Cartesian(double* Dist_Frac, double* Dist_Cart, int i){
	double x = Dist_Frac[0];
	double y = Dist_Frac[1];
	double z = Dist_Frac[2];
		
	x = x*a[i] + y*b[i]*temp_gamma_cos[i] + z*c[i]*temp_beta_cos[i];
	y = y*b[i]*temp_gamma_sin[i] + z*c[i]*temp_beta_term[i];
	z = z*c[i]*temp_gamma_term[i];
	Dist_Cart[0] = x; Dist_Cart[1] = y; Dist_Cart[2] = z;
}

// =============================================================================
// ----------------------------- Simple functions ------------------------------
// =============================================================================

// Return random number between 0 and 1
double Rand_0_1(){
	double x;
	do {
		x = rand();
		x = x/RAND_MAX;
	}
	while(x==0);
	
	return x;
}

// Calculates n!
double Facto(int n) {
	if (n < 2) 
		return 1.0;
  
  	else {
		double x = 1.0;
		for (int i = 2; i <= n; i++) {
			x *= i;
		}
		return x;
	}
}

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
		
		input >> n_frame >> n_mol;
		input >> snap_delay;
		input >> LAMBDA_I_H >> LAMBDA_I_E >> LAMBDA_S >> T >> H_OMEGA >> dist_tot >> n_try;
		//input >> n_mini_grid_a >> n_mini_grid_b >> n_mini_grid_c >> n_charges;
		input >> n_mini_grid_a >> n_mini_grid_b >> n_mini_grid_c >> tmp;
		input >> F_norm;
		
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
		cout << LAMBDA_I << " " << LAMBDA_S << " " << T << " " << H_OMEGA << " " << dist_tot << " " << n_try << " " << endl;
		cout << n_mini_grid_a << " " << n_mini_grid_b << " " << n_mini_grid_c << " " << n_charges << endl;
		cout << F_norm << endl;		

			for (int i=0; i<n_frame; i++){
				cout << "frame " << i << endl;
				
				for (int ii=0; ii<n_mol; ii++){
					input >> tmp >> mol_label_tmp >> n;
					mol_label.push_back(mol_label_tmp);
					cout << "molecule " << mol_label[ii] << " " << neigh_label[i][ii].size() << endl;
					
					for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
						cout << neigh_label[i][ii][jj] << " " << J_H[i][ii][jj] << " " << J_L[i][ii][jj] << endl;
					}
				}
			}
		}
		
	}
	else{
		cerr << "Error opening " << file_mc.str().c_str() << endl;
		exit(1);
	}
}

// Read .cell file
void Read_CELL(string input_file, string input_folder, bool print_results){
	double a_tmp, b_tmp, c_tmp, alpha_deg_tmp, beta_deg_tmp, gamma_deg_tmp, vol_box_tmp;
	double temp_alpha_cos_tmp, temp_beta_sin_tmp, temp_beta_cos_tmp, temp_gamma_sin_tmp, temp_gamma_cos_tmp, temp_beta_term_tmp, temp_gamma_term_tmp;

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
			input >> temp_alpha_cos_tmp >> temp_beta_sin_tmp >> temp_beta_cos_tmp >> temp_gamma_sin_tmp >> temp_gamma_cos_tmp >> temp_beta_term_tmp >> temp_gamma_term_tmp;
			
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
				cout << a[i] << " " << b[i] << " " << c[i] << " " << alpha_deg[i] << " " << beta_deg[i] << " " << gamma_deg[i] << endl;
				cout << temp_alpha_cos[i] << " " << temp_beta_sin[i] << " " << temp_beta_cos[i] << " " << temp_gamma_sin[i] << " " << temp_gamma_cos[i] << " " << temp_beta_term[i] << " " << temp_gamma_term[i] << endl;

			}
		}
		
	}
	else{
		cerr << "Error opening " << file_cell.str().c_str() << endl;
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
					cout << "molecule " << mol_label[ii] << " " << CM_x[i][ii] << " " << CM_y[i][ii] << " " << CM_z[i][ii] << endl;
				}
			}
		}
	}
	else{
		cerr << "Error opening " << file_cm.str().c_str() << endl;
		exit(1);
	}
}

// Read e_av part (average energies)
void Read_E_av(string input_file, string input_folder, bool print_results){
	double E_0_tmp, E_1_tmp;
	
	string tmp;
	stringstream file_e_av;
	file_e_av << input_folder << "/" << input_file << ".e_av";
	
	ifstream input(file_e_av.str().c_str(), ios::in);
	if (input){
		
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
	}
	else{
		cout << "[WARNING] Average energy file " << file_e_av.str().c_str() << " not found!\n[WARNING] E set to zero." << endl;
		
		for (int i=0; i<n_frame; i++) {
			E_0.push_back( vector<double> (n_mol, 0));
			E_1.push_back( vector<double> (n_mol, 0));
		}
	}
	
	// Print part
	if(print_results){
		for(int ii=0; ii<n_mol; ii++){
			cout << "molecule " << mol_label[ii] << " " << E_0[0][ii] << " " << E_1[0][ii] << endl;
		}
	}
}

// =============================================================================
// ------------------- Physical parameters between neighbors -------------------
// =============================================================================

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

// Calculates deltaE between molecules
void Calcul_DeltaE(bool print_results){
	cout << "[WARNING] The energies are supposed to be in kcal/mol, not eV!!!" << endl;
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
						//dE[i][ii].push_back(E_1[0][ii]+E_0[0][ll] - (E_0[0][ii]+E_1[0][ll]));
						dE[i][ii].push_back((E_1[0][ll]+E_0[0][ii] - (E_0[0][ll]+E_1[0][ii]))/23.06056);
						
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

double Calcul_DeltaV(int i, int mol_index_tmp, int neigh_index_tmp, int neigh_num_tmp, unsigned int charge_i_tmp, vector<int> curr_mol_tmp, vector<int> curr_grid_tmp){
		
	const double CST1 = 1.0/(4.0*PI*EPSILON_0*EPSILON_R);
	const double CUTOFF_2 = pow(CUTOFF_ELECTRO,2); //Cutoff electro square
	
	double dist, dist_neigh, dist_2, dist_neigh_2;
	double V_mol = 0.0, V_neigh = 0.0, V_electro = 0.0;
	
	double *CM_1_Cart, *CM_1_Frac, *CM_1_Neigh_Cart, *CM_1_Neigh_Frac, *CM_2_Cart, *CM_2_Frac, *Dist_Cart, *Dist_Frac;
	CM_1_Cart = new double[3];
	CM_1_Frac = new double[3];
	CM_1_Neigh_Cart = new double[3];
	CM_1_Neigh_Frac = new double[3];
	CM_2_Cart = new double[3];
	CM_2_Frac = new double[3];
	Dist_Cart = new double[3];
	Dist_Frac = new double[3];
	
	CM_1_Cart[0] = CM_x[i][mol_index_tmp];
	CM_1_Cart[1] = CM_y[i][mol_index_tmp];
	CM_1_Cart[2] = CM_z[i][mol_index_tmp];
	Cartesian_To_Fractional(CM_1_Cart, CM_1_Frac, i);
	
	CM_1_Frac[0] += double(box_a[curr_grid_tmp[charge_i_tmp]]);
	CM_1_Frac[1] += double(box_b[curr_grid_tmp[charge_i_tmp]]);
	CM_1_Frac[2] += double(box_c[curr_grid_tmp[charge_i_tmp]]);
	
	CM_1_Neigh_Cart[0] = CM_x[i][neigh_index_tmp];
	CM_1_Neigh_Cart[1] = CM_y[i][neigh_index_tmp];
	CM_1_Neigh_Cart[2] = CM_z[i][neigh_index_tmp];
	Cartesian_To_Fractional(CM_1_Neigh_Cart, CM_1_Neigh_Frac, i);
	
	CM_1_Neigh_Frac[0] += double(box_a[curr_grid_tmp[charge_i_tmp]]) + neigh_jump_vec_a[i][mol_index_tmp][neigh_num_tmp];
	CM_1_Neigh_Frac[1] += double(box_b[curr_grid_tmp[charge_i_tmp]]) + neigh_jump_vec_b[i][mol_index_tmp][neigh_num_tmp];
	CM_1_Neigh_Frac[2] += double(box_c[curr_grid_tmp[charge_i_tmp]]) + neigh_jump_vec_c[i][mol_index_tmp][neigh_num_tmp];
	
	for (unsigned int charge_i = 0; charge_i < curr_mol_tmp.size(); charge_i++){
		if (charge_i != charge_i_tmp){
			CM_2_Cart[0] = CM_x[i][curr_mol_tmp[charge_i]];
			CM_2_Cart[1] = CM_y[i][curr_mol_tmp[charge_i]];
			CM_2_Cart[2] = CM_z[i][curr_mol_tmp[charge_i]];
			Cartesian_To_Fractional(CM_2_Cart, CM_2_Frac, i);
			
			CM_2_Frac[0] += double(box_a[curr_grid_tmp[charge_i]]);
			CM_2_Frac[1] += double(box_b[curr_grid_tmp[charge_i]]);
			CM_2_Frac[2] += double(box_c[curr_grid_tmp[charge_i]]);
			
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
				
				V_electro += V_neigh - V_mol;
			}
		}
	}
	
	delete [] CM_1_Cart;
	delete [] CM_1_Frac;
	delete [] CM_1_Neigh_Cart;
	delete [] CM_1_Neigh_Frac;
	delete [] CM_2_Cart;
	delete [] CM_2_Frac;
	delete [] Dist_Cart;
	delete [] Dist_Frac;
	
	//cout << V_electro << endl;
	
	return V_electro;
}

// =============================================================================
// ------------------- Physical parameters between neighbors -------------------
// =============================================================================

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

double Marcus_Levich_Jortner_rate(double d_x_tmp, double d_y_tmp, double d_z_tmp, double dE_tmp, double J_H_tmp, double J_L_tmp){
	
	double dG0 = 0.0;	
	
	// CHECK SIGN!
	if (charge.compare("e") == 0)
		dG0 = -(d_x_tmp * F_x + d_y_tmp * F_y + d_z_tmp * F_z)*1E-8;
		
	else if (charge.compare("h") == 0)
		dG0 = (d_x_tmp * F_x + d_y_tmp * F_y + d_z_tmp * F_z)*1E-8;
		
	dG0 = dG0 + dE_tmp;
	
	int n = 0;
	double k_tmp = 0.0;
	double k_inter = MLJ_CST3[n]*exp(-pow(dG0 + LAMBDA_S + n*H_OMEGA,2)/(MLJ_CST2));
	
	while (k_inter > numeric_limits<double>::min()) {
		k_tmp += k_inter;
		n++;
		k_inter = MLJ_CST3[n]*exp(-pow(dG0 + LAMBDA_S + n*H_OMEGA,2)/(MLJ_CST2));
	}
	
	if (charge.compare("e") == 0)
		k_tmp = MLJ_CST1 * pow(J_L_tmp,2) * k_tmp;
	
	else if (charge.compare("h") == 0)
		k_tmp = MLJ_CST1 * pow(J_H_tmp,2) * k_tmp;
		
	return k_tmp;
	
}

double Marcus_Levich_Jortner_rate_electro(int i, int mol_index_tmp, int neigh_index_tmp, int neigh_num_tmp, double d_x_tmp, double d_y_tmp, double d_z_tmp, double dE_tmp, double J_H_tmp, double J_L_tmp, vector<int> curr_mol_tmp, vector<int> curr_grid_tmp, unsigned int charge_i_tmp){
	
	double dV = 0.0;
	double dG0 = 0.0;	
	
	// Calcul DeltaV
	dV = Calcul_DeltaV(i, mol_index_tmp, neigh_index_tmp, neigh_num_tmp, charge_i_tmp, curr_mol_tmp, curr_grid_tmp);

	// CHECK SIGN!
	if (charge.compare("e") == 0)
		dG0 = -(d_x_tmp * F_x + d_y_tmp * F_y + d_z_tmp * F_z)*1E-8;
		
	else if (charge.compare("h") == 0)
		dG0 = (d_x_tmp * F_x + d_y_tmp * F_y + d_z_tmp * F_z)*1E-8;
	
	//printf("%e %e\n", dG0, dV);	
	
	dG0 = dG0 + dE_tmp + dV;
	
	int n = 0;
	double k_tmp = 0.0;
	double k_inter = MLJ_CST3[n]*exp(-pow(dG0 + LAMBDA_S + n*H_OMEGA,2)/(MLJ_CST2));
	
	while (k_inter > numeric_limits<double>::min()) {
		k_tmp += k_inter;
		n++;
		k_inter = MLJ_CST3[n]*exp(-pow(dG0 + LAMBDA_S + n*H_OMEGA,2)/(MLJ_CST2));
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
				
				double k_tmp = Marcus_Levich_Jortner_rate(d_x[i][ii][jj], d_y[i][ii][jj], d_z[i][ii][jj], dE[i][ii][jj], J_H[i][ii][jj], J_L[i][ii][jj]);
				
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
						
						neigh_jump_vec_a[i][ii].insert(neigh_jump_vec_a[i][ii].begin(), -neigh_jump_vec_a[i][ll][jj]);
						neigh_jump_vec_b[i][ii].insert(neigh_jump_vec_b[i][ii].begin(), -neigh_jump_vec_b[i][ll][jj]);
						neigh_jump_vec_c[i][ii].insert(neigh_jump_vec_c[i][ii].begin(), -neigh_jump_vec_c[i][ll][jj]);

						dE[i][ii].insert(dE[i][ii].begin(), -dE[i][ll][jj]);
						
						J_H[i][ii].insert(J_H[i][ii].begin(), J_H[i][ll][jj]);
						J_L[i][ii].insert(J_L[i][ii].begin(), J_L[i][ll][jj]);
							
						double k_tmp = Marcus_Levich_Jortner_rate(d_x[i][ii].front(), d_y[i][ii].front(), d_z[i][ii].front(), dE[i][ii].front(), J_H[i][ii].front(), J_L[i][ii].front());
						
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

// =============================================================================
// ------------------------------- Grid building -------------------------------
// =============================================================================

void Build_Grid(bool print_results) {
	
	int n_grid = n_mini_grid_a*n_mini_grid_b*n_mini_grid_c;
	
	for (int a=0; a<n_mini_grid_a; a++){
		for (int b=0; b<n_mini_grid_b; b++){
			for (int c=0; c<n_mini_grid_c; c++){
				box_a.push_back(a);
				box_b.push_back(b);
				box_c.push_back(c);

			}
		}
	}
	
	for (int x=0; x<n_grid; x++){
		box_neigh_label.push_back( vector< int > () );
		box_neigh_a.push_back( vector< int > () );
		box_neigh_b.push_back( vector< int > () );
		box_neigh_c.push_back( vector< int > () );
		
		for (int y=0; y<n_grid; y++){
			if ((abs(box_a[y]-box_a[x]) < 2 && abs(box_b[y]-box_b[x]) < 2 && abs(box_c[y]-box_c[x]) < 2) && (abs(box_a[y]-box_a[x]) != 0 || abs(box_b[y]-box_b[x]) != 0 || abs(box_c[y]-box_c[x]) != 0)){
				box_neigh_a[x].push_back( box_a[y]-box_a[x] );
				box_neigh_b[x].push_back( box_b[y]-box_b[x] );
				box_neigh_c[x].push_back( box_c[y]-box_c[x] );
				box_neigh_label[x].push_back( y );
				
			}
		}
	}
	
	if (print_results){
		
		for (unsigned int x=0; x<box_neigh_label.size(); x++){
			cout << "mini-grid " << x << endl;			
			for (unsigned int xx=0; xx<box_neigh_label[x].size(); xx++){
				cout << box_neigh_label[x][xx] << " " << box_neigh_a[x][xx] << " " << box_neigh_b[x][xx] << " " << box_neigh_c[x][xx] << endl;
			}
		}
	}
}


// =============================================================================
// --------------------------------- Print part --------------------------------
// =============================================================================

// Print all the informations used to calculate k
void Print_Summary(string output_folder) {
	stringstream OUT_TOT, T, P;
	double uF_x, uF_y, uF_z;
  
	//T << theta_deg;
	//P << phi_deg;
	//OUT_TOT << output_folder << "/info_tot_" << T.str().c_str() << "_" << P.str().c_str() << ".dat";

	OUT_TOT << output_folder << "/info_tot_" << charge.c_str() << "_" << F_dir.c_str() << ".dat";

	uF_x = F_x/F_norm;
	uF_y = F_y/F_norm;
	uF_z = F_z/F_norm;

	FILE * pFile;

	pFile = fopen(OUT_TOT.str().c_str(), "w");
	if (pFile==NULL) {
		int wait = 0; 
		while (wait<10 && pFile==NULL){
			cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) << " seconds to write a file" << endl;
			usleep(10*(wait+1)*(wait+1));
			pFile = fopen(OUT_TOT.str().c_str(), "w");
			wait++;
		}
		if (wait==10 && pFile==NULL){
			cerr << "[ERROR] Impossible to write " << OUT_TOT.str().c_str() << "! Exiting..." << endl;
			exit (1);
		}
	}
	fprintf(pFile,"Electric field unit vectors: (%f, %f, %f)\n\n", uF_x, uF_y, uF_z);
	for (int i=0; i<n_frame; i++){
		fprintf(pFile,"Frame %d\n", i);
		for (int ii=0; ii<n_mol; ii++){
			fprintf(pFile,"Molecule %d | %d neighbors\n", mol_label[ii], int(neigh_label[i][ii].size()));
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
				fprintf(pFile,"%6d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %14.5e\n", neigh_label[i][ii][jj], d_x[i][ii][jj], d_y[i][ii][jj], d_z[i][ii][jj], dE[i][ii][jj], J_H[i][ii][jj], J_L[i][ii][jj], k[i][ii][jj]);
			}
		}
	}

	fclose(pFile);
}

// =============================================================================
// ----------------------- Monte-Carlo related functions -----------------------
// =============================================================================

// Choose a starting molecule randomly
int Choose_Mol_RND(int frame){ 
	int mol; double k_sum;
	
	do {
		mol = rand()%n_mol;
		k_sum = 0.0;
		
		for(unsigned int jj=0; jj<neigh_label[frame][mol].size(); jj++) {
			k_sum = k_sum + k_inv[frame][mol][jj];
		}
	}
	while(neigh_label[frame][mol].size()==0 || neigh_label[frame][mol].size() == 0);

	return mol;
}

void Dispatch_Mol_RND(int frame, vector< vector<bool> > grid_occ, int *pos){ 
	
	int pos_a=0, pos_b=0, pos_c=0, mol=0, box=0; 
	
	//double k_sum;
	
	do {
		pos_a = rand()%n_mini_grid_a;
		pos_b = rand()%n_mini_grid_b;
		pos_c = rand()%n_mini_grid_c;

		mol = rand()%n_mol;
		
		//k_sum = 0.0;
		//
		//for(unsigned int jj=0; jj<neigh_label[frame][mol].size(); jj++) {
		//	k_sum = k_sum + k_inv[frame][mol][jj];
		//}
		
		// Find the number of the corresponding mini-grid
		for (unsigned int x=0; x<box_neigh_label.size(); x++){
			if (pos_a == box_a[x] && pos_b == box_b[x] && pos_c == box_c[x]){
				box = x;
				break;
			}
		}
	}
	while(neigh_label[frame][mol].size() == 0 || grid_occ[mol][box] == true);
	
	pos[0] = mol;
	pos[1] = box;
}

void Dispatch_Mol_begin(int frame, vector< vector<bool> > grid_occ, int *pos){ 
	
	int pos_a=0, pos_b=0, pos_c=0, mol=0, box=0; 
	
	//double k_sum;
	
	do {
		if (F_dir.compare("a") == 0) {
			pos_a = 0;
			pos_b = rand()%n_mini_grid_b;
			pos_c = rand()%n_mini_grid_c;
		}
		
		else if (F_dir.compare("b") == 0) {
			pos_a = rand()%n_mini_grid_a;
			pos_b = 0;
			pos_c = rand()%n_mini_grid_c;
		}
		
		else if (F_dir.compare("c") == 0) {
			pos_a = rand()%n_mini_grid_a;
			pos_b = rand()%n_mini_grid_b;
			pos_c = 0;
		}	
		
		mol = rand()%n_mol;
		
		//k_sum = 0.0;
		//
		//for(unsigned int jj=0; jj<neigh_label[frame][mol].size(); jj++) {
		//	k_sum = k_sum + k_inv[frame][mol][jj];
		//}
		
		// Find the number of the corresponding mini-grid
		for (unsigned int x=0; x<box_neigh_label.size(); x++){
			if (pos_a == box_a[x] && pos_b == box_b[x] && pos_c == box_c[x]){
				box = x;
				break;
			}
		}
	}
	while(neigh_label[frame][mol].size() == 0 || grid_occ[mol][box] == true);
	
	pos[0] = mol;
	pos[1] = box;
}

// BKL algorithm
void MC_BKL(string output_folder){
	
	cout << "[INFO] Using the BKL algorithm." << endl;
	
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
	
	// Start the BKL algorithm
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
		fprintf(pFile,"Frame = %d\n", i);
		fprintf(pFile,"Electric Field Unit Vectors: (%f, %f, %f)\n", uF_x, uF_y, uF_z);
		fprintf(pFile,"Number of Charges = %d\n", n_charges);
		fprintf(pFile,"Density of Charges = %.5e charges/cm3\n", double(n_charges)/(vol_box[i]*n_mini_grid_a*n_mini_grid_b*n_mini_grid_c*10e-24));
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
			Dispatch_Mol_RND(i, grid_occ, pos);
			
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
						
						double k_tmp = numeric_limits<double>::min();
						if (tmp_neigh_index != tmp_mol_index){
							k_tmp = Marcus_Levich_Jortner_rate_electro(i, tmp_mol_index, tmp_neigh_index, tmp_neigh_num, d_x[i][tmp_mol_index][jj], d_y[i][tmp_mol_index][jj], d_z[i][tmp_mol_index][jj], dE[i][tmp_mol_index][jj], J_H[i][tmp_mol_index][jj], J_L[i][tmp_mol_index][jj], curr_mol, curr_grid, charge_i);
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
			double random_k = Rand_0_1() * sum_k;
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
						//printf("out of system\n");
					}
					
					// Do not jump if the next molecule is occupied
					else if (grid_occ[tmp_curr_mol][tmp_curr_grid] == true){
						previous_jump_ok = false;
						exit_loop = true;
						//printf("occupied\n");
					}
					
					else{
						
						// Calculate the total time
						total_time_try += -log(Rand_0_1())/(sum_k);
						
						// Calculate the distance traveled by the charge and the total distance
						double event_dist = (d_x[i][event_mol_index[event]][event_neigh_num[event]] * uF_x + d_y[i][event_mol_index[event]][event_neigh_num[event]] * uF_y + d_z[i][event_mol_index[event]][event_neigh_num[event]] * uF_z)*1E-8;
						dist[event_charge[event]] += event_dist;
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
						if ((F_dir.compare("a") == 0 && tmp_curr_grid_a >= n_mini_grid_a-1) || (F_dir.compare("b") == 0 && tmp_curr_grid_b >= n_mini_grid_b-1) || (F_dir.compare("c") == 0 && tmp_curr_grid_c >= n_mini_grid_c-1)){
							
							// The charge is removed
							// curr_mol.erase(curr_mol.begin()+event_charge[event]);
							// curr_grid.erase(curr_grid.begin()+event_charge[event]);
							
							// The charge appears at the beginning of the system
							int *pos;
							pos = new int[2];
							
							Dispatch_Mol_begin(i, grid_occ, pos);
							
							curr_grid[event_charge[event]] = pos[1];
							curr_mol[event_charge[event]] = pos[0];
							grid_occ[pos[0]][pos[1]] = true;
							
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
					}
				}
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
		fprintf(pFile,"Number of tries = %d\n", n_try);
		fprintf(pFile,"Average Time = %e\n", total_time_try/double(n_try));
		fprintf(pFile,"Average Distance = %e\n", total_dist_try/double(n_try));
		fprintf(pFile,"Mobility of the Frame %d = %lf\n", i, mu_frame.back());
		fclose(pFile);
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
	fprintf(pFile,"\n==================\n==================\n");
	fprintf(pFile,"mu_av = %lf\n", mu_moy);
	fclose(pFile);
	
}

// BKL algorithm
void MC_BKL_MT(string output_folder){
	
	cout << "[INFO] Using the BKL algorithm with multithreading." << endl;
	
	// Average mobility for each frame
	vector<double> mu_frame (n_frame, 0.0);
	
	#pragma omp parallel for private(grid_occ)	
	
	// Start the BKL algorithm
	for (int i=0; i<n_frame; i++){
		
		// Variables for each charge
		vector<int> curr_mol, curr_grid; // Number of the molecule in the mini-grid, and number of this mini-grid
		vector<double> dist, jump; // Distance and number of jumps for each charge
		
		// Variables for a chosen event
		vector<double> event_k; // Transfer rate
		vector<int> event_charge, event_mol_index, event_neigh_num, event_neigh_index; // Initial molecule and neighbor
		
		// Total time and distance
		double total_time_try, total_dist_try;

		double uF_x, uF_y, uF_z;
		uF_x = F_x/F_norm;
		uF_y = F_y/F_norm;
		uF_z = F_z/F_norm;
		
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
		fprintf(pFile,"Frame = %d\n", i);
		fprintf(pFile,"Electric Field Unit Vectors: (%f, %f, %f)\n", uF_x, uF_y, uF_z);
		fprintf(pFile,"Number of Charges = %d\n", n_charges);
		fprintf(pFile,"Density of Charges = %.5e charges/cm3\n", double(n_charges)/(vol_box[i]*n_mini_grid_a*n_mini_grid_b*n_mini_grid_c*10e-24));
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
			Dispatch_Mol_RND(i, grid_occ, pos);
			
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
						
						double k_tmp = numeric_limits<double>::min();
						if (tmp_neigh_index != tmp_mol_index){
							k_tmp = Marcus_Levich_Jortner_rate_electro(i, tmp_mol_index, tmp_neigh_index, tmp_neigh_num, d_x[i][tmp_mol_index][jj], d_y[i][tmp_mol_index][jj], d_z[i][tmp_mol_index][jj], dE[i][tmp_mol_index][jj], J_H[i][tmp_mol_index][jj], J_L[i][tmp_mol_index][jj], curr_mol, curr_grid, charge_i);
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
			double random_k = Rand_0_1() * sum_k;
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
						//printf("out of system\n");
					}
					
					// Do not jump if the next molecule is occupied
					else if (grid_occ[tmp_curr_mol][tmp_curr_grid] == true){
						previous_jump_ok = false;
						exit_loop = true;
						//printf("occupied\n");
					}
					
					else{
						
						// Calculate the total time
						total_time_try += -log(Rand_0_1())/(sum_k);
						
						// Calculate the distance traveled by the charge and the total distance
						double event_dist = (d_x[i][event_mol_index[event]][event_neigh_num[event]] * uF_x + d_y[i][event_mol_index[event]][event_neigh_num[event]] * uF_y + d_z[i][event_mol_index[event]][event_neigh_num[event]] * uF_z)*1E-8;
						dist[event_charge[event]] += event_dist;
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
						if ((F_dir.compare("a") == 0 && tmp_curr_grid_a >= n_mini_grid_a-1) || (F_dir.compare("b") == 0 && tmp_curr_grid_b >= n_mini_grid_b-1) || (F_dir.compare("c") == 0 && tmp_curr_grid_c >= n_mini_grid_c-1)){
							
							// The charge is removed
							// curr_mol.erase(curr_mol.begin()+event_charge[event]);
							// curr_grid.erase(curr_grid.begin()+event_charge[event]);
							
							// The charge appears at the beginning of the system
							int *pos;
							pos = new int[2];
							
							Dispatch_Mol_begin(i, grid_occ, pos);
							
							curr_grid[event_charge[event]] = pos[1];
							curr_mol[event_charge[event]] = pos[0];
							grid_occ[pos[0]][pos[1]] = true;
							
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
					}
				}
			}
		}

		// Calculates the final mobility for the frame
		mu_frame[i] = total_dist_try/(total_time_try*F_norm);

		// Writes a summary for the frame
		pFile=fopen(OUT_SIMU.str().c_str(), "w");
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
		fprintf(pFile,"Number of tries = %d\n", n_try);
		fprintf(pFile,"Average Time = %e\n", total_time_try/double(n_try));
		fprintf(pFile,"Average Distance = %e\n", total_dist_try/double(n_try));
		fprintf(pFile,"Mobility of the Frame %d = %lf\n", i, mu_frame[i]);
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
	fprintf(pFile,"mu_av = %lf\n", mu_moy);
	fclose(pFile);
	
}

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
		fprintf(pFile,"Frame = %d\n", i);
		fprintf(pFile,"Electric Field Unit Vectors: (%f, %f, %f)\n", uF_x, uF_y, uF_z);
		fprintf(pFile,"Number of Charges = %d\n", n_charges);
		fprintf(pFile,"Density of Charges = %.5e charges/cm3\n", double(n_charges)/(vol_box[i]*n_mini_grid_a*n_mini_grid_b*n_mini_grid_c*10e-24));
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
			Dispatch_Mol_RND(i, grid_occ, pos);
			
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
				if ((F_dir.compare("a") == 0 && tmp_curr_grid_a >= n_mini_grid_a-1) || (F_dir.compare("b") == 0 && tmp_curr_grid_b >= n_mini_grid_b-1) || (F_dir.compare("c") == 0 && tmp_curr_grid_c >= n_mini_grid_c-1)){
					
					// The charge is removed
					// curr_mol.erase(curr_mol.begin()+event_charge[event]);
					// curr_grid.erase(curr_grid.begin()+event_charge[event]);
					
					// The charge appears at the beginning of the system
					int *pos;
					pos = new int[2];
					
					Dispatch_Mol_begin(i, grid_occ, pos);
					
					curr_grid[event_charge[event]] = pos[1];
					curr_mol[event_charge[event]] = pos[0];
					grid_occ[pos[0]][pos[1]] = true;
					
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
		fprintf(pFile,"Number of tries = %d\n", n_try);
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
	fprintf(pFile,"\n==================\n==================\n");
	fprintf(pFile,"mu_av = %lf\n", mu_moy);
	fclose(pFile);
	
}

// FRM algorithm
void MC_FRM_MT(string output_folder){
	
	cout << "[INFO] Using the FRM algorithm with multithread." << endl;
	
	// Average mobility for each frame
	vector<double> mu_frame (n_frame, 0.0);
	
	#pragma omp parallel for private(grid_occ)
	
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

		double uF_x, uF_y, uF_z;
		uF_x = F_x/F_norm;
		uF_y = F_y/F_norm;
		uF_z = F_z/F_norm;
		
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
		fprintf(pFile,"Frame = %d\n", i);
		fprintf(pFile,"Electric Field Unit Vectors: (%f, %f, %f)\n", uF_x, uF_y, uF_z);
		fprintf(pFile,"Number of Charges = %d\n", n_charges);
		fprintf(pFile,"Density of Charges = %.5e charges/cm3\n", double(n_charges)/(vol_box[i]*n_mini_grid_a*n_mini_grid_b*n_mini_grid_c*10e-24));
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
			Dispatch_Mol_RND(i, grid_occ, pos);
			
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
				if ((F_dir.compare("a") == 0 && tmp_curr_grid_a >= n_mini_grid_a-1) || (F_dir.compare("b") == 0 && tmp_curr_grid_b >= n_mini_grid_b-1) || (F_dir.compare("c") == 0 && tmp_curr_grid_c >= n_mini_grid_c-1)){
					
					// The charge is removed
					// curr_mol.erase(curr_mol.begin()+event_charge[event]);
					// curr_grid.erase(curr_grid.begin()+event_charge[event]);
					
					// The charge appears at the beginning of the system
					int *pos;
					pos = new int[2];
					
					Dispatch_Mol_begin(i, grid_occ, pos);
					
					curr_grid[event_charge[event]] = pos[1];
					curr_mol[event_charge[event]] = pos[0];
					grid_occ[pos[0]][pos[1]] = true;
					
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
		pFile=fopen(OUT_SIMU.str().c_str(), "w");
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
		fprintf(pFile,"Number of tries = %d\n", n_try);
		fprintf(pFile,"Average Time = %e\n", total_time_try/double(n_try));
		fprintf(pFile,"Average Distance = %e\n", total_dist_try/double(n_try));
		fprintf(pFile,"Mobility of the Frame %d = %lf\n", i, mu_frame[i]);
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
	fprintf(pFile,"mu_av = %lf\n", mu_moy);
	fclose(pFile);
	
}

// =============================================================================
// --------------------------------- Clear part --------------------------------
// =============================================================================

void Clear_All(){
	a.clear(); b.clear(); c.clear(); alpha_deg.clear(); beta_deg.clear(); gamma_deg.clear(); vol_box.clear(); 
	temp_alpha_cos.clear(); temp_beta_sin.clear(); temp_beta_cos.clear(); temp_gamma_sin.clear(); temp_gamma_cos.clear(); temp_beta_term.clear(); temp_gamma_term.clear(); 

	mol_label.clear();
	CM_x.clear(); CM_y.clear(); CM_z.clear(); 
	E_0.clear(); E_1.clear();

	box_a.clear(); box_b.clear(); box_c.clear(); 
	grid_occ.clear();
	box_neigh_a.clear(); box_neigh_b.clear(); box_neigh_c.clear(); box_neigh_label.clear(); 

	neigh_label.clear();
	d_x.clear(); d_y.clear(); d_z.clear();
	dE.clear();
	J_H.clear(); J_L.clear();
	neigh_jump_vec_a.clear(); neigh_jump_vec_b.clear(); neigh_jump_vec_c.clear(); 
	k.clear(); k_inv.clear();

	MLJ_CST3.clear();

	event_k_1.clear(); event_k_2.clear(); 
	event_charge_1.clear(); event_mol_index_1.clear(); event_neigh_num_1.clear(); event_neigh_index_1.clear(); event_charge_2.clear(); event_mol_index_2.clear(); event_neigh_num_2.clear(), event_neigh_index_2.clear();
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

	int s;
	string input_file;
	string input_folder = ".";
	string output_folder = ".";
	string method;
	
  	while ((s = getopt_long (argc, argv, "I:i:o:c:d:n:m:", NULL, NULL)) != -1){
      	switch (s){
			case 'I':
				input_file = optarg;
				break;
				
			case 'i':
				input_folder = optarg;
				struct stat st_in;
				if(stat(input_folder.c_str(), &st_in) != 0){
					cout << "[ERROR] Input folder " << input_folder << " doesn't exist! Exiting..." << endl;
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
			}
	}
	
	// Check that the charge and the direction of the electric field are specified
	if (charge.compare("e") == 0 || charge.compare("h") == 0)
		cout << "[INFO] The charge is: " << charge << endl;
		
	else {
		cerr << "[ERROR] Charge not specified! Please use -c {e,h}. Exiting..." << endl;
		exit(1);
	}
	
	if (F_dir.compare("a") == 0 || F_dir.compare("b") == 0 || F_dir.compare("c") == 0)
		cout << "[INFO] Electric field is along the '" << F_dir << "' direction." << endl;
		
	else {
		cerr << "[ERROR] Direction of the electric field not specified! Please use -d {a,b,c}. Exiting..." << endl;
		exit(1);
	}

	// Read the required files
	Read_MC(input_file, input_folder, false);
	Read_CELL(input_file, input_folder, false);
	Read_CM(input_file, input_folder, false);
	Read_E_av(input_file, input_folder, false);
	
	// Calculates distances and DeltaE
	Calcul_Dist(false);
	Calcul_DeltaE(false);

	// Build the grid
	Build_Grid(false);
	
	// Save the triangular matrix
	vector< vector< vector<int> > > neigh_label_ref = neigh_label;
	vector< vector< vector<int> > > neigh_jump_vec_a_ref = neigh_jump_vec_a, neigh_jump_vec_b_ref = neigh_jump_vec_b, neigh_jump_vec_c_ref = neigh_jump_vec_c;
	vector< vector< vector<double> > > J_H_ref = J_H, J_L_ref = J_L;
	vector< vector< vector<double> > > d_x_ref = d_x, d_y_ref = d_y, d_z_ref = d_z;
	vector< vector< vector<double> > > dE_ref = dE;				
		
//	for (int i=90; i<91; i=i+15){
//		for (int j=0; j<1; j=j+15){
			
//			cout << "[INFO] Running simulation for phi = " << j << endl;

			// theta_deg = i; phi_deg = j;
			
			// theta_rad = 2 * PI * i; theta_rad = theta_rad/360.0;
      		// phi_rad = 2 * PI * j; phi_rad = phi_rad/360.0;
      		
			// F_x = sin(theta_rad); F_x = F_x * cos(phi_rad); 
			// F_y = sin(theta_rad); F_y = F_y * sin(phi_rad); 
			// F_z = cos(theta_rad); 

			// Calculates the electric field unit vector
			double *F_tmp_frac, *F_tmp_cart;
			F_tmp_frac = new double[3];
			F_tmp_cart = new double[3];

			if (F_dir.compare("a") == 0) {
				if (charge.compare("e") == 0) {
					F_tmp_frac[0] = 1.0;
					LAMBDA_I = LAMBDA_I_E;
				}
				else {
					F_tmp_frac[0] = -1.0;
					LAMBDA_I = LAMBDA_I_H;
				}
				F_tmp_frac[1] = 0.0;
				F_tmp_frac[2] = 0.0;
			}
				
			else if (F_dir.compare("b") == 0) {
				F_tmp_frac[0] = 0.0;
				if (charge.compare("e") == 0) {
					F_tmp_frac[1] = 1.0;
					LAMBDA_I = LAMBDA_I_E;
				}
				else {
					F_tmp_frac[1] = -1.0;
					LAMBDA_I = LAMBDA_I_H;
				}
				F_tmp_frac[2] = 0.0;
			}
				
			else if (F_dir.compare("c") == 0) {
				F_tmp_frac[0] = 0.0;
				F_tmp_frac[1] = 0.0;
				if (charge.compare("e") == 0) {
					F_tmp_frac[2] = 1.0;
					LAMBDA_I = LAMBDA_I_E;
				}
				else {
					F_tmp_frac[2] = -1.0;
					LAMBDA_I = LAMBDA_I_H;
				}
			}
						
			Fractional_To_Cartesian(F_tmp_frac, F_tmp_cart, 0);
			
			double F_norm_tmp = sqrt(pow(F_tmp_cart[0],2) + pow(F_tmp_cart[1],2) + pow(F_tmp_cart[2],2));
			F_x = F_tmp_cart[0]/F_norm_tmp;
			F_y = F_tmp_cart[1]/F_norm_tmp;
			F_z = F_tmp_cart[2]/F_norm_tmp;
			
			delete [] F_tmp_frac; delete [] F_tmp_cart;
			
			if (fabs(F_x) < 1E-10) F_x = 0.0;
			if (fabs(F_y) < 1E-10) F_y = 0.0;
			if (fabs(F_z) < 1E-10) F_z = 0.0;
			
			F_x = F_x * F_norm; 
			F_y = F_y * F_norm; 
			F_z = F_z * F_norm;
			
			// Calculate transfer rates and the full matrix, mostly for information
			Marcus_Levich_Jortner_CST(); // Calculates constants
			Calcul_k(false);
			Full_Matrix();
	
			// Print a summary
			Print_Summary(output_folder);
			
			// Clear some tables
			// J_H.clear();
			// J_L.clear();
			// dE.clear();
			
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
			
//		}
//	}

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
