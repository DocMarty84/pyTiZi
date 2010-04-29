/**
 *******************************************************************************
 * Copyright (C) 2010 Nicolas Martinelli, nicolas.martinelli@gmail.com         *
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

#include <iostream> //Entrées-sorties standard
#include <fstream> //Entrées-sorties sur fichiers
#include <string> //Chaines de caracteres
#include <limits> //Pour aller à la fin d'une ligne, par exemple
#include <iomanip> //Manipulation des sorties, par exemple format scientifique
#include <sstream> //Pour les conversion entre types de variables
#include <vector>
#include <math.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <algorithm>

using namespace std;

double K_BOLTZ = 8.617343e-5; //Boltzmann constant in eV
double H_BAR = 6.58211899e-16; // h/2PI in eV
double PI = 3.141592653589793;

vector< vector<double> > CM_x, CM_y, CM_z; // Center of masses

bool pbc[3]; // Periodic boundary conditions
vector<double> a, b, c, alpha_deg, beta_deg, gamma_deg; // Cell parameters
vector<double> temp_alpha_cos, temp_beta_sin, temp_beta_cos, temp_gamma_sin, temp_gamma_cos, temp_beta_term, temp_gamma_term; // Parameters for fractional coordinates

int n_frame, n_mol; // Number of frame and molecule
double snap_delay;
double LAMBDA_I, LAMBDA_S, T, H_OMEGA, dist_tot;
unsigned int n_try, n_charges;
double F_norm;
vector<int> mol_label;
vector< vector< vector<int> > > neigh_label;
vector< vector< vector<double> > > J_H;
vector< vector< vector<double> > > J_L;

vector< vector<double> > E_0, E_1;

vector< vector< vector<double> > > d_x, d_y, d_z;
vector< vector< vector<double> > > dE;
vector< vector< vector<double> > > k, k_inv;
string charge;

double theta_deg, phi_deg, theta_rad, phi_rad;
double F_x, F_y, F_z;

int n_mini_grid_a;
int n_mini_grid_b;
int n_mini_grid_c;

// Vector for the change in mini-grid
vector< vector< vector<int> > > neigh_jump_vec_a;
vector< vector< vector<int> > > neigh_jump_vec_b;
vector< vector< vector<int> > > neigh_jump_vec_c;

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
		input >> LAMBDA_I >> LAMBDA_S >> T >> H_OMEGA >> dist_tot >> n_try;
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
	double a_tmp, b_tmp, c_tmp, alpha_deg_tmp, beta_deg_tmp, gamma_deg_tmp;
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
			input >> a_tmp >> b_tmp >> c_tmp >> alpha_deg_tmp >> beta_deg_tmp >> gamma_deg_tmp;
			input >> temp_alpha_cos_tmp >> temp_beta_sin_tmp >> temp_beta_cos_tmp >> temp_gamma_sin_tmp >> temp_gamma_cos_tmp >> temp_beta_term_tmp >> temp_gamma_term_tmp;
			
			a.push_back(a_tmp); 
			b.push_back(b_tmp); 
			c.push_back(c_tmp);
			alpha_deg.push_back(alpha_deg_tmp); 
			beta_deg.push_back(beta_deg_tmp); 
			gamma_deg.push_back(gamma_deg_tmp);

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
			input >> tmp >> tmp >> E_0_tmp >> E_1_tmp;

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
	
	vector< vector< vector<int> > > neigh_jump_vec_a;
	vector< vector< vector<int> > > neigh_jump_vec_b;
	vector< vector< vector<int> > > neigh_jump_vec_c;
	
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
						dE[i][ii].push_back(E_1[0][ll]+E_0[0][ii] - (E_0[0][ll]+E_1[0][ii]));
						
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

// =============================================================================
// ------------------- Physical parameters between neighbors -------------------
// =============================================================================

double Marcus_Levich_Jortner_rate(double d_x_tmp, double d_y_tmp, double d_z_tmp, double dE_tmp, double J_H_tmp, double J_L_tmp){
	
	double S, CST1, CST2;
	S = LAMBDA_I/H_OMEGA;
	CST1 = (2*PI/H_BAR)*(1.0/(sqrt(4*PI*LAMBDA_S*K_BOLTZ*T)));
	CST2 = 4*LAMBDA_S*K_BOLTZ*T;
	
	double dG0 = 0.0;	
	double k_tmp = 0.0;
	
	// CHECK SIGN!
	if (charge.compare("e") == 0)
		dG0 = -(d_x_tmp * F_x + d_y_tmp * F_y + d_z_tmp * F_z)*1E-8;
		
	if (charge.compare("h") == 0)
		dG0 = (d_x_tmp * F_x + d_y_tmp * F_y + d_z_tmp * F_z)*1E-8;
		
	dG0 = dG0 + dE_tmp;
	
	for (int n=0; n<50; n++){
		k_tmp = k_tmp + exp(-S)*(pow(S,n)/Facto(n))*exp(-pow(dG0 + LAMBDA_S + n*H_OMEGA,2)/(CST2));
	}
	
	if (charge.compare("e") == 0)
		k_tmp = CST1 * pow(J_L_tmp,2) * k_tmp;
	
	if (charge.compare("h") == 0)
		k_tmp = CST1 * pow(J_H_tmp,2) * k_tmp;
		
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
// --------------------------------- Print part --------------------------------
// =============================================================================

// Print all the informations used to calculate k
void Print_Summary(string output_folder) {
	stringstream OUT_TOT, T, P;
	double uF_x, uF_y, uF_z;
  
	T << theta_deg;
	P << phi_deg;
	OUT_TOT << output_folder << "/info_tot_" << T.str().c_str() << "_" << P.str().c_str() << ".dat";

	uF_x = F_x/F_norm;
	uF_y = F_y/F_norm;
	uF_z = F_z/F_norm;

	FILE * pFile;

	pFile = fopen(OUT_TOT.str().c_str(), "w");
	if (pFile==NULL) { 
		cerr << "[ERROR] Impossible to write " << OUT_TOT.str().c_str() << "! Exiting..." << endl;
		exit (1);
	}

	fprintf(pFile,"Electric field unit vectors: (%f, %f, %f)\n\n", uF_x, uF_y, uF_z);
	for (int i=0; i<n_frame; i++){
		fprintf(pFile,"Frame %d\n", i);
		for (int ii=0; ii<n_mol; ii++){
			fprintf(pFile,"Molecule %d | %d neighbors\n", mol_label[ii], neigh_label[i][ii].size());
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

void Dispatch_Mol_RND(int frame, vector< vector< vector< vector<bool> > > > grid_occ, int *pos){ 
	
	int pos_a, pos_b, pos_c, mol; 
	
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
	}
	while(neigh_label[frame][mol].size() == 0 || grid_occ[pos_a][pos_b][pos_c][mol] == true);
	
	pos[0] = pos_a;
	pos[1] = pos_b;
	pos[2] = pos_c;
	pos[3] = mol;
	
}

// BKL algorithm
void MC_BKL(string output_folder){
	
	// Grid of charges
	vector< vector< vector< vector<bool> > > > grid_occ;
	
	// Variables for each charge
	vector<int> pos_a, pos_b, pos_c; // Position in the grid
	vector<int> curr_mol; // Number of the molecule in the mini-grid
	vector<double> dist, jump; // Distance and number of jumps for each charge
	
	// Variables for a chosen event
	vector<double> event_k; // Transfer rate
	vector<int> event_mol_init, event_mol_fin; // Initial molecule and neighbor
	
	// Total time and distance, for each try
	vector<double> total_time_try, total_dist_try;
	
	// Average mobility for each frame
	vector<double> mu_frame;
	
	double uF_x, uF_y, uF_z;
	uF_x = F_x/F_norm;
	uF_y = F_y/F_norm;
	uF_z = F_z/F_norm;
	
	// Check if it's possible to write output files
	stringstream OUT_SIMU, OUT_ERROR, T, P;
  
	T << theta_deg;
	P << phi_deg;
	OUT_SIMU << output_folder << "/simu_" << T.str().c_str() << "_" << P.str().c_str() << ".out";
	OUT_ERROR << output_folder << "/error_" << T.str().c_str() << "_" << P.str().c_str() << ".out";
  
	FILE * pFile;

	pFile = fopen(OUT_SIMU.str().c_str(), "w");
	if (pFile==NULL) { 
		cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
		exit (1);
	}
	fclose(pFile);
	
	// Start the BKL algorithm
	for (int i=0; i<n_frame; i++){
		
		// Clear the total distance and time of the previous frame
		total_time_try.clear();
		total_dist_try.clear();
				
		for (unsigned int charge_try=0; charge_try<n_try; charge_try++){
			
			// Generate the grid
			grid_occ.clear();
					
			for (int a=0; a<n_mini_grid_a; a++){
				grid_occ.push_back( vector< vector< vector<bool> > > ());
				
				for (int b=0; b<n_mini_grid_b; b++){
					grid_occ[a].push_back( vector< vector<bool> > ());
				
					for (int c=0; c<n_mini_grid_c; c++){
						grid_occ[a][b].push_back( vector<bool> ());
				
						for (int n=0; n<n_mol; n++){
							grid_occ[a][b][c].push_back( false );
				
						}
					}
				}
			}

			// Put the charges in the grid
			pos_a.clear(); pos_b.clear(); pos_c.clear();
			curr_mol.clear();

			int *pos;
			pos = new int[4];
			for (unsigned int charge_i; charge_i<n_charges; charge_i++){
				Dispatch_Mol_RND(i, grid_occ, pos);
				
				pos_a.push_back(pos[0]);
				pos_b.push_back(pos[1]);
				pos_c.push_back(pos[2]);
				curr_mol.push_back(pos[3]);
			}
			delete [] pos;

			// Set distances, number of jumps and the travel time to zero
			total_time_try.push_back(0.0);
			total_dist_try.push_back(0.0);
			
			dist.clear();
			jump.clear(); 
			for (unsigned int charge_i; charge_i<n_charges; charge_i++){
				dist.push_back(0.0);
				jump.push_back(0.0);
			}
			
			bool previous_jump_ok = true;
			bool exit_loop = false;
			double sum_k = 0.0;
			
			while (curr_mol.size()>0){
				
				// Calculates the transfer rates
				if (previous_jump_ok){
					
					// List all the possible events
					event_k.clear();
					event_mol_init.clear();
					event_mol_fin.clear();
					
					for (unsigned int charge_i; charge_i<curr_mol.size(); charge_i++){
						for (unsigned int jj=0; jj<neigh_label[i][curr_mol[charge_i]].size(); jj++){
							double k_tmp = Marcus_Levich_Jortner_rate(d_x[i][curr_mol[charge_i]][jj], d_y[i][curr_mol[charge_i]][jj], d_z[i][curr_mol[charge_i]][jj], dE[i][curr_mol[charge_i]][jj], J_H[i][curr_mol[charge_i]][jj], J_L[i][curr_mol[charge_i]][jj], pos_a[charge_i], pos_b[charge_i], pos_c[charge_i]);
							
							event_k.push_back(k_tmp);
							event_mol_init.push_back(charge_i);
							event_mol_fin.push_back(jj);
							
						}
					}
					
					// Sum of transfer rates of all the events
					sum_k = 0.0;
					for (unsigned int event=0; event<event_k.size(); event++){
						sum_k += event_k[event];
					}
				}
				
				else{
					exit_loop = false;
				}
				
				// Choose an event randomly
				double random_k = Rand_0_1() * sum_k;
				double partial_sum_k = 0.0;
				
				// Find this event in the list
				exit_loop = false;
				
				for (unsigned int event=0; event<event_k.size() && exit_loop == false; event++){
					partial_sum_k += event_k[event];
					
					if (partial_sum_k > random_k){
						
						// Calculate the new position in the grid
						int tmp_pos_a = pos_a[event_mol_init[event]] + neigh_jump_vec_a[i][curr_mol[event_mol_init[event]]][event_mol_fin[event]];
						int tmp_pos_b = pos_b[event_mol_init[event]] + neigh_jump_vec_b[i][curr_mol[event_mol_init[event]]][event_mol_fin[event]];
						int tmp_pos_c = pos_c[event_mol_init[event]] + neigh_jump_vec_c[i][curr_mol[event_mol_init[event]]][event_mol_fin[event]];
						
						// No PBC in this case
						//tmp_pos_a = tmp_pos_a % n_mini_grid_a;
						//tmp_pos_b = tmp_pos_b % n_mini_grid_b;
						//tmp_pos_c = tmp_pos_c % n_mini_grid_c;
						
						// Search the corresponding next molecule
						int tmp_curr_mol = 0;
						for (int ii=0; ii<n_mol; ii++){
							if (neigh_label[i][curr_mol[event_mol_init[event]]][event_mol_fin[event]] == mol_label[ii]){
								tmp_curr_mol = ii;
								break;
							}
						}
						
						// Do not jump if the charge goes out of the system
						if (tmp_pos_a < 0 || tmp_pos_b < 0 || tmp_pos_c < 0){
							previous_jump_ok = false;
							exit_loop = true;
						}
						
						// Do not jump if the next molecule is occupied
						else if (grid_occ[abs(tmp_pos_a)][abs(tmp_pos_b)][abs(tmp_pos_c)][tmp_curr_mol] == true){
							previous_jump_ok = false;
							exit_loop = true;
						}
						
						else{
							
							// Calculate the total time
							total_time_try.back() += -(sum_k)*log(Rand_0_1());
							
							// Calculate the distance traveled by the charge and the total distance
							double event_dist = (d_x[i][curr_mol[event_mol_init[event]]][event_mol_fin[event]] * uF_x + d_y[i][curr_mol[event_mol_init[event]]][event_mol_fin[event]] * uF_y + d_z[i][curr_mol[event_mol_init[event]]][event_mol_fin[event]] * uF_z)*1E-8;
							dist[event_mol_init[event]] += event_dist;
							total_dist_try.back() += event_dist;
							
							// Calculate the number of jumps for each charge
							jump[event_mol_init[event]] += 1.0;
							
							// Set the occupancy of the grid
							grid_occ[pos_a[event_mol_init[event]]][pos_b[event_mol_init[event]]][pos_c[event_mol_init[event]]][curr_mol[event_mol_init[event]]] = false;
							
							// Set the new position in the grid from temp values
							pos_a[event_mol_init[event]] = tmp_pos_a;
							pos_b[event_mol_init[event]] = tmp_pos_b;
							pos_c[event_mol_init[event]] = tmp_pos_c;
							curr_mol[event_mol_init[event]] = tmp_curr_mol;
							
							// Check if the charge arrived at the end of the grid
							if (pos_a[event_mol_init[event]] >= n_mini_grid_a){
								
								// The charge is removed
								pos_a.erase(pos_a.begin()+event_mol_init[event]);
								pos_b.erase(pos_b.begin()+event_mol_init[event]);
								pos_c.erase(pos_c.begin()+event_mol_init[event]);
								curr_mol.erase(curr_mol.begin()+event_mol_init[event]);
								
							}
							else{
								grid_occ[pos_a[event_mol_init[event]]][pos_b[event_mol_init[event]]][pos_c[event_mol_init[event]]][curr_mol[event_mol_init[event]]] = true;
							}
							
							previous_jump_ok = true;
							exit_loop = true;
						}
					}
				}
			}
			
			// Print summary for the current try
			pFile=fopen(OUT_SIMU.str().c_str(), "a");
			if (pFile==NULL) { 
				cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
				exit (1);
			}
			fprintf(pFile,"Frame = %d\n", i);
			fprintf(pFile,"Try = %d\n", charge_try);
			fprintf(pFile,"Electric Field Unit Vectors: (%f, %f, %f)\n", uF_x, uF_y, uF_z);
			fprintf(pFile,"Number of Charges = %d\n", n_charges);
			fprintf(pFile,"Total Time = %e\n", total_time_try.back());
			fprintf(pFile,"Total Distance = %e\n", total_dist_try.back());
			fprintf(pFile,"Mu_try = %lf\n", total_dist_try.back()/(total_time_try.back()*F_norm));
			fclose(pFile);			

		}

		// Calculates the average travel time and distance for the frame
		double total_time_av = 0.0; 
		double total_dist_av = 0.0; 
		double dbl_n_try = n_try;
		
		for(unsigned int charge_try=0; charge_try<n_try; charge_try++) {
			total_time_av = total_time_av + total_time_try[charge_try];
			total_dist_av = total_dist_av + total_dist_try[charge_try];
		}
		total_time_av = total_time_av/dbl_n_try;
		total_dist_av = total_dist_av/dbl_n_try;
		
		mu_frame.push_back(total_dist_av/(total_time_av*F_norm));

		// Writes a summary for the frame
		pFile=fopen(OUT_SIMU.str().c_str(), "a");
		if (pFile==NULL) { 
			cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
			exit (1);
		}
		fprintf(pFile,"Frame = %d\n", i);
		fprintf(pFile,"Electric field unit vectors: (%f, %f, %f)\n", uF_x, uF_y, uF_z);
		fprintf(pFile,"Number of tries = %d\n", n_try);
		fprintf(pFile,"Average Time = %e\n", total_time_av);
		fprintf(pFile,"Average Distance = %e\n", total_dist_av);
		fprintf(pFile,"Mobility of the Frame = %lf\n", mu_frame.back());
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
		cerr << "[ERROR] Impossible to write " << OUT_SIMU.str().c_str() << "! Exiting..." << endl;
		exit (1);
	}
	fprintf(pFile,"\n==================\n==================\n");
	fprintf(pFile,"mu_av = %lf\n", mu_moy);
	fclose(pFile);
	
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
	
  	while ((s = getopt_long (argc, argv, "I:i:o:c:", NULL, NULL)) != -1){
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
			}
	}
	
	// Check that the charge is specified
	if (charge.compare("e") == 0)
		cout << "[INFO] The charge is an electron." << endl;
		
	else if (charge.compare("h") == 0)
		cout << "[INFO] The charge is a hole." << endl;
		
	else {
		cerr << "[ERROR] Charge (-c e or -c h) not specified! Exiting..." << endl;
		exit(1);
	}
	
	// Read the required files
	Read_MC(input_file, input_folder, false);
	Read_CELL(input_file, input_folder, false);
	Read_CM(input_file, input_folder, false);
	Read_E_av(input_file, input_folder, false);
	
	// Calcul distances and DeltaE
	Calcul_Dist(false);
	Calcul_DeltaE(false);
	
	// Save the triangular matrix
	vector< vector< vector<int> > > neigh_label_ref = neigh_label;
	vector< vector< vector<double> > > J_H_ref = J_H, J_L_ref = J_L;
	vector< vector< vector<double> > > d_x_ref = d_x, d_y_ref = d_y, d_z_ref = d_z;
	vector< vector< vector<double> > > dE_ref = dE;
	
	for (int i=90; i<91; i=i+15){
		for (int j=0; j<360; j=j+15){
			
			cout << "[INFO] Running simulation for phi = " << j << endl;
			
			theta_deg = i; phi_deg = j;
			
			theta_rad = 2 * PI * i; theta_rad = theta_rad/360.0;
      		phi_rad = 2 * PI * j; phi_rad = phi_rad/360.0;
      		
			F_x = sin(theta_rad); F_x = F_x * cos(phi_rad); 
			F_y = sin(theta_rad); F_y = F_y * sin(phi_rad); 
			F_z = cos(theta_rad); 
			
			if (fabs(F_x) < 1E-10) F_x = 0.0;
			if (fabs(F_y) < 1E-10) F_y = 0.0;
			if (fabs(F_z) < 1E-10) F_z = 0.0;
			
			F_x = F_x * F_norm; 
			F_y = F_y * F_norm; 
			F_z = F_z * F_norm;
			
			// Calculate transfer rates and the full matrix
			Calcul_k(false);
			Full_Matrix();
			
			// Print a summary
			Print_Summary(output_folder);
			
			// Clear some tables
			J_H.clear();
			J_L.clear();
			dE.clear();
			
			// Calculates 1/k and clear the k table
			Inverse_Clear_k(false);
			k.clear();
			
			MC_BKL(output_folder);
			
			// Clear everything (not necessary) and set to reference values
			neigh_label.clear(); neigh_label = neigh_label_ref;
			J_H.clear(); J_H = J_H_ref;
			J_L.clear(); J_L = J_L_ref;
			d_x.clear(); d_x = d_x_ref;
			d_y.clear(); d_y = d_y_ref;
			d_z.clear(); d_z = d_z_ref;
			dE.clear();	dE = dE_ref;
			k.clear();
			k_inv.clear();
			
		}
	}
	
	// Get stop time
	t_stop = time(NULL);
	timeinfo_stop = asctime(localtime(&t_stop));

	// Print final information
	cout << "[INFO] Start time: " << timeinfo_start;	
	cout << "[INFO] Stop time: " << timeinfo_stop;
	cout << "[INFO] Calculation took " << t_stop-t_start << " seconds." << endl << endl;
	cout << "===============================================================================" << endl;
	cout << "-------------------- Exiting normally, everything was ok! ---------------------" << endl;
	cout << "===============================================================================" << endl;
	
	return 0;
}