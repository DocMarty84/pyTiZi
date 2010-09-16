// C++ libraries
#include <iostream> 
#include <fstream> 
#include <string> 
#include <limits> 
#include <iomanip> 
#include <sstream>
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

double ***x_cart, ***y_cart, ***z_cart; // Cartesian coordinates 
int ***atomic_number; // Atomic numbers
double ***atomic_mass; // Atomic masses
int ***atomic_valence; // Atomic valences
char ***symbol; // Atomic symbols
double **CM_x, **CM_y, **CM_z; // Center of masses
int *n_electrons; // Number of valence electrons in a molecule
int n_frame, n_mol, *n_atom, *mol_label; // Number of frame, molecule, atom list and label of molecules
bool *J; // Calculate neighbors of this variable?

bool pbc[3]; // Periodic boundary conditions
double cutoff; // Cutoff distance
double *a, *b, *c, *alpha_deg, *beta_deg, *gamma_deg; // Cell parameters
double *temp_alpha_cos, *temp_beta_sin, *temp_beta_cos, *temp_gamma_sin, *temp_gamma_cos, *temp_beta_term, *temp_gamma_term; // Parameters for fractional coordinates

vector< vector< vector<int> > > neigh_label;
vector< vector< vector<int> > > neigh_jump_vec_a, neigh_jump_vec_b, neigh_jump_vec_c; // Vector for the change in mini-grid

bool sign;
int coeff_H_lign, coeff_H_row, coeff_L_lign, coeff_L_row; 

vector< vector< vector<double> > > J_H, J_L;

bool MT = 1;

void Read_XYZ(string input_file, string input_folder, bool print_results){
	string tmp;
	stringstream file_xyz;
	file_xyz << input_folder << "/" << input_file << ".xyz";
	
	ifstream input(file_xyz.str().c_str(), ios::in);
	if (input){
		input >> n_frame >> n_mol;
		
		n_atom = new int[n_mol];
		for(int ii=0; ii<n_mol; ii++){
			input >> n_atom[ii];
		}
		
		x_cart = new double**[n_frame];
		y_cart = new double**[n_frame];
		z_cart = new double**[n_frame];
		atomic_number = new int**[n_frame];
		atomic_mass = new double**[n_frame];
		atomic_valence = new int**[n_frame];
		symbol = new char**[n_frame];
		mol_label = new int[n_mol];
		
		for(int i=0; i<n_frame; i++){
			x_cart[i] = new double*[n_mol];
			y_cart[i] = new double*[n_mol];
			z_cart[i] = new double*[n_mol];
			atomic_number[i] = new int*[n_mol];
			atomic_mass[i] = new double*[n_mol];
			atomic_valence[i] = new int*[n_mol];
			symbol[i] = new char*[n_mol];
			
			for(int ii=0; ii<n_mol; ii++){
				x_cart[i][ii] = new double[n_atom[ii]];
				y_cart[i][ii] = new double[n_atom[ii]];
				z_cart[i][ii] = new double[n_atom[ii]];
				atomic_number[i][ii] = new int[n_atom[ii]];
				atomic_mass[i][ii] = new double[n_atom[ii]];
				atomic_valence[i][ii] = new int[n_atom[ii]];
				symbol[i][ii] = new char[n_atom[ii]];
			}
		}
		
		for(int i=0; i<n_frame; i++){
			input >> tmp >> tmp;
			
			for(int ii=0; ii<n_mol; ii++){
				input >> tmp >> mol_label[ii];
				
				for(int iii=0; iii<n_atom[ii]; iii++){
					input >> symbol[i][ii][iii] >> atomic_number[i][ii][iii] >> atomic_mass[i][ii][iii] >> atomic_valence[i][ii][iii] >> x_cart[i][ii][iii] >> y_cart[i][ii][iii] >> z_cart[i][ii][iii];
				}
			}
		}
		
		if (print_results){
			for(int i=0; i<n_frame; i++){
				cout << "frame " << i << endl;
				
				for(int ii=0; ii<n_mol; ii++){
					cout << "molecule " << mol_label[ii] << endl;
					
					for(int iii=0; iii<n_atom[ii]; iii++){
						cout << symbol[i][ii][iii] << " " << atomic_number[i][ii][iii] << " " << atomic_mass[i][ii][iii] << " " << atomic_valence[i][ii][iii] << " " << x_cart[i][ii][iii] << " " << y_cart[i][ii][iii] << " " << z_cart[i][ii][iii] << endl;;
					}
				}
			}
		}
		
		input.close();
	}
	
	else{
		cerr << "Error opening " << file_xyz.str().c_str() << endl;
		exit(1);
	}
}

void Read_CELL(string input_file, string input_folder, bool print_results){
	string tmp;
	stringstream file_cell;
	file_cell << input_folder << "/" << input_file << ".cell";
	
	ifstream input(file_cell.str().c_str(), ios::in);
	if (input){
		a = new double[n_frame];
		b = new double[n_frame];
		c = new double[n_frame];
		alpha_deg = new double[n_frame];
		beta_deg = new double[n_frame];
		gamma_deg = new double[n_frame];
		
		temp_alpha_cos = new double[n_frame]; 
		temp_beta_sin = new double[n_frame];
		temp_beta_cos = new double[n_frame];
		temp_gamma_sin = new double[n_frame];
		temp_gamma_cos = new double[n_frame];
		temp_beta_term = new double[n_frame];
		temp_gamma_term = new double[n_frame];

		input >> tmp >> pbc[0] >> pbc[1] >> pbc[2];
		input >> tmp >> cutoff;
		for(int i=0; i<n_frame; i++){
			input >> tmp >> tmp;
			input >> a[i] >> b[i] >> c[i] >> alpha_deg[i] >> beta_deg[i] >> gamma_deg[i];
			input >> temp_alpha_cos[i] >> temp_beta_sin[i] >> temp_beta_cos[i] >> temp_gamma_sin[i] >> temp_gamma_cos[i] >> temp_beta_term[i] >> temp_gamma_term[i];
		}
		
		input.close();
		
		if (print_results){
			cout << "PBC " << pbc[0] << " " << pbc[1] << " " << pbc[2] << endl;
			cout << "cutoff " << cutoff << endl;
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

void Read_CM(string input_file, string input_folder, bool print_results){
	string tmp;
	stringstream file_cm;
	file_cm << input_folder << "/" << input_file << ".cm";
	
	ifstream input(file_cm.str().c_str(), ios::in);
	if (input){
		
		CM_x = new double*[n_frame];
		CM_y = new double*[n_frame];
		CM_z = new double*[n_frame];
		n_electrons = new int[n_mol];
		mol_label = new int[n_mol];
		J = new bool[n_mol];
		
		for(int i=0; i<n_frame; i++){
			CM_x[i] = new double[n_mol];
			CM_y[i] = new double[n_mol];
			CM_z[i] = new double[n_mol];
		}
		
		for(int i=0; i<n_frame; i++){
			input >> tmp >> tmp;
			for(int ii=0; ii<n_mol; ii++){
				input >> tmp >> mol_label[ii] >> n_electrons[ii] >> CM_x[i][ii] >> CM_y[i][ii] >> CM_z[i][ii] >> J[ii];
			}
		}
		
		input.close();
		
		if(print_results){
			for(int i=0; i<n_frame; i++){
				cout << "frame " << i << endl;
				for(int ii=0; ii<n_mol; ii++){
					cout << "molecule " << mol_label[ii] << " " << n_electrons[ii] << " " << CM_x[i][ii] << " " << CM_y[i][ii] << " " << CM_z[i][ii] << " " << J[ii] << endl;
				}
			}
		}
	}
	else{
		cerr << "Error opening " << file_cm.str().c_str() << endl;
		exit(1);
	}
}

void Read_ZIN(string input_file, string input_folder, bool print_results){
	string tmp;
	stringstream file_zin;
	file_zin << input_folder << "/" << input_file << ".zin";
	
	ifstream input(file_zin.str().c_str(), ios::in);
	if (input){
		input >> sign >> coeff_H_lign >> coeff_H_row >> coeff_L_lign >> coeff_L_row; 
		input.close();
		
		if(print_results){
			cout << "Read sign: " << sign << endl;
			cout << "coeff_H_lign: " << coeff_H_lign << endl;
			cout << "coeff_H_row: " << coeff_H_row << endl;
			cout << "coeff_L_lign: " << coeff_L_lign << endl;
			cout << "coeff_L_row: " << coeff_L_row << endl;
		}
	}
	else{
		cerr << "Error opening " << file_zin.str().c_str() << endl;
		exit(1);
	}
}

void Read_NB(string input_file, string input_folder, bool print_results){
	string tmp;
	stringstream file_nb;
	file_nb << input_folder << "/" << input_file << ".nb";
	
	//unsigned int tmp_n_neigh;
	//int tmp_neigh_label;
	
	ifstream input(file_nb.str().c_str(), ios::in);
	if (input){
		
		neigh_label.clear();
		input >> n_frame >> n_mol;

		/*
		for (int i=0; i<n_frame; i++){
			
			neigh_label.push_back( vector< vector<int> > ());
			input >> tmp >> tmp;
			
			for (int ii=0; ii<n_mol; ii++){
				
				neigh_label[i].push_back( vector<int> ());
				input >> tmp >> tmp >> tmp_n_neigh;
				
				for (unsigned int jj=0; jj<tmp_n_neigh; jj++){
					input >> tmp_neigh_label;
					neigh_label[i][ii].push_back(tmp_neigh_label);
				}
			}
		}
		*/
		input.close();
		
		if (print_results){
			cout << n_frame << " " << n_mol << endl;
			for (int i=0; i<n_frame; i++){
				cout << "frame " << i << endl;
				for (int ii=0; ii<n_mol; ii++){
					cout << "molecule " << ii << " " << neigh_label[i][ii].size() << endl;
					for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
						cout << neigh_label[i][ii][jj] << endl;
					}
				}
			}
		}
	}
		
	else{
		cerr << "Error opening " << file_nb.str().c_str() << endl;
		exit(1);
	}
}

void Read_FULL(string input_file, string result_folder){
	string tmp;
	stringstream file_full;
	file_full << result_folder << "/" << input_file << ".full";
	
	ifstream input(file_full.str().c_str(), ios::in);
	
	if (input){
		int f, mol1, mol2, mol1_index;
		double j_h, j_l, h_1, h_2, l_1, l_2;
		
		neigh_label.clear();
		J_H.clear(); J_L.clear(); 
		
		for(int i=0; i<n_frame; i++){
			neigh_label.push_back( vector< vector<int> > ());
			J_H.push_back( vector< vector<double> > ());
			J_L.push_back( vector< vector<double> > ());
			
			for(int ii=0; ii<n_mol; ii++){
				neigh_label[i].push_back( vector<int> ());
				J_H[i].push_back( vector<double> ());
				J_L[i].push_back( vector<double> ());
			}
		}
		
		while (!input.eof()){
			input >> f >> mol1 >> mol2 >> j_h >> h_1 >> h_2 >> j_l >> l_1 >> l_2;
			for (int ii=0; ii<n_mol; ii++){
				if (mol_label[ii]==mol1){
					mol1_index = ii;
					break;
				}
			}
			
			// Fill the vector and sort it at the same time
			if (neigh_label[f][mol1_index].size() == 0){
				neigh_label[f][mol1_index].push_back(mol2);
				J_H[f][mol1_index].push_back(1000*j_h*(h_1*h_2/(fabs(h_1*h_2))));
				J_L[f][mol1_index].push_back(1000*j_l*(l_1*l_2/(fabs(l_1*l_2))));
			}

			else if (mol2 > neigh_label[f][mol1_index].back()) {
				neigh_label[f][mol1_index].push_back(mol2);
				J_H[f][mol1_index].push_back(1000*j_h*(h_1*h_2/(fabs(h_1*h_2))));
				J_L[f][mol1_index].push_back(1000*j_l*(l_1*l_2/(fabs(l_1*l_2))));
			}
			
			else if (mol2 < neigh_label[f][mol1_index].front()) {
				neigh_label[f][mol1_index].insert(neigh_label[f][mol1_index].begin(), mol2);
				J_H[f][mol1_index].insert(J_H[f][mol1_index].begin(), 1000*j_h*(h_1*h_2/(fabs(h_1*h_2))));
				J_L[f][mol1_index].insert(J_L[f][mol1_index].begin(), 1000*j_l*(l_1*l_2/(fabs(l_1*l_2))));
			}
			
			else {
				vector<int>::iterator it_n = neigh_label[f][mol1_index].begin();
				vector<double>::iterator it_h = J_H[f][mol1_index].begin();
				vector<double>::iterator it_l = J_L[f][mol1_index].begin();
				for (it_n=neigh_label[f][mol1_index].begin(); it_n<neigh_label[f][mol1_index].end(); it_n++){
					if (mol2<*(it_n+1) && mol2>*(it_n)) {
						neigh_label[f][mol1_index].insert(it_n+1, mol2);
						J_H[f][mol1_index].insert(it_h+1, (1000*j_h*(h_1*h_2/(fabs(h_1*h_2)))));
						J_L[f][mol1_index].insert(it_l+1, (1000*j_l*(l_1*l_2/(fabs(l_1*l_2)))));
						break;
					}
					it_h++;
					it_l++;
				}
			}
		}
		input.close();
	}
	
	else{
		cerr << "Error opening " << file_full.str().c_str() << endl;
		exit(1);			
	}
}

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

void Mol_Frac_To_Cart(double** Coord_Frac, double** Coord_Cart, int frame, int atom){
	for (int i=0; i<atom; i++){
		Fractional_To_Cartesian(Coord_Frac[i], Coord_Cart[i], frame);
	}
}

void Mol_Cart_To_Frac(double** Coord_Cart, double** Coord_Frac, int frame, int atom){
	for (int i=0; i<atom; i++){
		Cartesian_To_Fractional(Coord_Cart[i], Coord_Frac[i], frame);
	}
}

void Find_Neighbors_Sphere_MT(string input_file, string output_folder, bool print_results){

	neigh_label.clear();
	neigh_jump_vec_a.clear(); neigh_jump_vec_b.clear(); neigh_jump_vec_c.clear();
	
	for(int i=0; i<n_frame; i++){
		neigh_label.push_back( vector< vector<int> > ());
		neigh_jump_vec_a.push_back( vector< vector<int> > ());
		neigh_jump_vec_b.push_back( vector< vector<int> > ());
		neigh_jump_vec_c.push_back( vector< vector<int> > ());
		
		for(int ii=0; ii<n_mol; ii++){
			neigh_label[i].push_back( vector<int> ());
			neigh_jump_vec_a[i].push_back( vector<int> ());
			neigh_jump_vec_b[i].push_back( vector<int> ());
			neigh_jump_vec_c[i].push_back( vector<int> ());
		}
	}

	double CutOff_square = pow(cutoff, 2);
	
	#pragma omp parallel for
	
	// Start calculation
	for (int i=0; i<n_frame; i++){
		
		double Dist_Norm_square = 0.0;
		double *Dist_Cart, *Dist_Frac;
		int *tmp_neigh_jump_vec;
		Dist_Cart = new double[3];
		Dist_Frac = new double[3];
		tmp_neigh_jump_vec = new int[3];
		
		for (int ii=0; ii<n_mol; ii++){
			if(J[ii]){
				for (int jj=0; jj<n_mol-(ii+1); jj++){
					Dist_Cart[0] = CM_x[i][jj+(ii+1)] - CM_x[i][ii];
					Dist_Cart[1] = CM_y[i][jj+(ii+1)] - CM_y[i][ii];
					Dist_Cart[2] = CM_z[i][jj+(ii+1)] - CM_z[i][ii];
					
					Cartesian_To_Fractional(Dist_Cart, Dist_Frac, i);

					// Periodic boundary conditions calculations
					for (int k=0; k<3; k++){
						tmp_neigh_jump_vec[k] = 0;
						if (fabs(Dist_Frac[k]) > 0.5 && pbc[k]){
							if (Dist_Frac[k] < 0.0){
								Dist_Frac[k] = Dist_Frac[k] + 1.0;
								tmp_neigh_jump_vec[k] = 1;
							}
							else{
								Dist_Frac[k] = Dist_Frac[k] - 1.0;
								tmp_neigh_jump_vec[k] = -1;
							}
						}
					}
					
					Fractional_To_Cartesian(Dist_Frac, Dist_Cart, i);
					
					Dist_Norm_square = pow(Dist_Cart[0],2) + pow(Dist_Cart[1],2) + pow(Dist_Cart[2],2);
					
					if (Dist_Norm_square < CutOff_square && Dist_Norm_square != 0){
						neigh_label[i][ii].push_back(mol_label[jj+(ii+1)]);
						neigh_jump_vec_a[i][ii].push_back(tmp_neigh_jump_vec[0]);
						neigh_jump_vec_b[i][ii].push_back(tmp_neigh_jump_vec[1]);
						neigh_jump_vec_c[i][ii].push_back(tmp_neigh_jump_vec[2]);
						
						stringstream output_filename;
						stringstream s_frame, s_mol_n1, s_mol_n2;
						
						//s_frame << i;
						//s_mol_n1 << mol_label[ii];
						//s_mol_n2 << neigh_label[i][ii].back();
						//output_filename << output_folder.c_str() << "/frame_" << s_frame.str().c_str() << "/dimer_" << s_mol_n1.str().c_str() << "_" << s_mol_n2.str().c_str() << ".dist";
						
						//ofstream output(output_filename.str().c_str(), ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier
						//if (output){
						//	output << Dist_Cart[0] << " " << Dist_Cart[1] << " " << Dist_Cart[2] <<endl;
						//}
						//else
						//	cerr << "Error opening " << output_filename.str().c_str() << endl;	
					}
				}
			}
		}
		
		delete [] Dist_Cart; delete [] Dist_Frac;
	}
	
	#pragma omp barrier
	
	if(print_results){
		for (int i=0; i<n_frame; i++){
			cout << "frame " << i << endl;
			for (int ii=0; ii<n_mol; ii++){
				cout << "Molecule " << mol_label[ii] << ": " << neigh_label[i][ii].size() << " Neighbors\n";
				for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
					cout << "Vector " << neigh_jump_vec_a[i][ii][jj] << " " << neigh_jump_vec_b[i][ii][jj] << " " << neigh_jump_vec_c[i][ii][jj] << endl;
				}
			}
		}	
	}
}

void Write_DAT(string input_file, double **mol1, double **mol2, int frame, int mol_n1, int mol_n2){
	
	stringstream output_filename;
	stringstream s_frame, s_mol_n1, s_mol_n2;
	
	s_frame << frame;
	s_mol_n1 << mol_label[mol_n1];
	s_mol_n2 << mol_label[mol_n2];
	output_filename << input_file.c_str() << "/frame_" << s_frame.str().c_str() << "/cell_" << s_mol_n1.str().c_str() << "_" << s_mol_n2.str().c_str() << ".dat";
	
	FILE * pFile;

	pFile = fopen (output_filename.str().c_str(),"w");
	
	fprintf(pFile, "%d %d\n", n_atom[mol_n1], n_atom[mol_n2]);
	for (int iii=0; iii<n_atom[mol_n1]; iii++)
		fprintf(pFile, "%10.5f%10.5f%10.5f%4d\n", mol1[iii][0], mol1[iii][1], mol1[iii][2], atomic_number[frame][mol_n1][iii]);
	
	for (int jjj=0; jjj<n_atom[mol_n2]; jjj++)
		fprintf(pFile, "%10.5f%10.5f%10.5f%4d\n", mol2[jjj][0], mol2[jjj][1], mol2[jjj][2], atomic_number[frame][mol_n2][jjj]);
	
	fclose (pFile);

	
}

void Write_INP(string input_file, double **mol, int frame, int mol_n){
	
	int n_s=0, n_p=0, n_d=0;
	int n_exc=50;
	
	for (int iii=0; iii<n_atom[mol_n]; iii++){
		if (atomic_number[frame][mol_n][iii] == 1){
			n_s = n_s++;
		}
		
		else if (atomic_number[frame][mol_n][iii] == 44){
			n_d = n_d++;
		}
		
		else{
			n_p = n_p++;
		}
	}

	stringstream output_filename;
	stringstream s_frame, s_mol_n;
	
	s_frame << frame;
	s_mol_n << mol_label[mol_n];
	output_filename << input_file.c_str() << "/frame_" << s_frame.str().c_str() << "/molecule_" << s_mol_n.str().c_str() << ".inp";
		
	FILE * pFile;

	pFile = fopen (output_filename.str().c_str(),"w");

	fprintf(pFile, " $TITLEI\n");
	fprintf(pFile, "INPUT CREATED BY Mon cul c'est du poulet\n\n");
	fprintf(pFile, " $END\n\n");
	fprintf(pFile, " $CONTRL\n\n");
	fprintf(pFile, " SCFTYP     RHF      RUNTYP     ENERGY     ENTTYP     COORD    UNITS      ANGS\n\n");
	fprintf(pFile, "INTTYP          1   IAPX            3   NAT          %6d   NEL       %4d \n", n_atom[mol_n], n_electrons[mol_n]);
	fprintf(pFile, "MULT            1   ITMAX          50 \n\n");
	fprintf(pFile, " SCFTOL   0.0000010  \n\n");
	fprintf(pFile, "! ***** Basis set and C.I. size inFormation *****\n\n");
	fprintf(pFile, "   DYNAL(1) =    0%5d%5d%5d  0  4000%5d \n\n", n_s, n_p, n_d, n_electrons[mol_n]/2);
	fprintf(pFile, "! ***** Interaction factors *****\n\n");
	fprintf(pFile, "   INTFA(1) =   1.00000  1.26700  0.58500  1.00000  1.00000  1.00000\n\n");
	fprintf(pFile, "! ***** Output file name *****\n\n");
	fprintf(pFile, "   ONAME = file99\n\n");
	fprintf(pFile, " $END\n\n");
	fprintf(pFile, " $DATAIN\n");
	for (int iii=0; iii<n_atom[mol_n]; iii++){
		fprintf(pFile, "%10.5f%10.5f%10.5f%4d\n", mol[iii][0], mol[iii][1], mol[iii][2], atomic_number[frame][mol_n][iii]);
	}
	fprintf(pFile, "\n");
	fprintf(pFile, " $END\n");
	fprintf(pFile, " $CIINPU\n\n");
	fprintf(pFile, "! ***** Configuration Interaction specification *****\n\n");
	fprintf(pFile, "    2    1 %2d    1    0    0    0    1    1    1 %2d\n", n_exc, n_exc); // A faire
	fprintf(pFile, "  -80000.0  0.000000\n\n");
	fprintf(pFile, "00001%5d%5d\n", n_electrons[mol_n]/2, n_electrons[mol_n]/2);
	fprintf(pFile, "00021%5d%5d%5d%5d\n\n\n\n", (n_electrons[mol_n]/2)+1, n_electrons[mol_n]/2, (n_electrons[mol_n]/2)+1, n_electrons[mol_n]/2);
	fprintf(pFile, " $END");
		
	fclose (pFile);
}

void Write_morange(string input_file, int frame){
	
	stringstream output_filename;
	stringstream s_frame;
	
	s_frame << frame;
	output_filename << input_file.c_str() << "/frame_" << s_frame.str().c_str() << "/morange.inp";
	ofstream output(output_filename.str().c_str(), ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier
	
	if (output){
		output << n_electrons[0]/2 - 1 << " " << n_electrons[0]/2 + 1 << endl;
		output << n_electrons[0]/2 - 1 << " " << n_electrons[0]/2 + 1;
		output.close();
	}
	else 
		cerr << "Error opening " << output_filename.str().c_str() << endl;
}

void Write_CMD(string input_file, string zindo_folder, string output_folder, string log_file, int frame, int mol_n1, int mol_n2){

	stringstream output_filename;
	stringstream s_frame, s_mol_n1, s_mol_n2;
	
	s_frame << frame;
	s_mol_n1 << mol_label[mol_n1];
	s_mol_n2 << mol_label[mol_n2];
	output_filename << input_file.c_str() << "/frame_" << s_frame.str().c_str() << "/run_" << s_mol_n1.str().c_str() << "_" << s_mol_n2.str().c_str() << ".cmd";
	ofstream output(output_filename.str().c_str(), ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier
	
	if (output){
		output << "#!/bin/bash" << endl << endl;
		output << "cp cell_"  <<  mol_label[mol_n1] << "_" << mol_label[mol_n2] << ".dat cell.dat" << endl << endl;
		output << "nlaya=1" << endl;
		output << "nlayb=1" << endl;
		output << "nlayc=1" << endl << endl;
		output << "a=1" << endl;
		output << "b=1" << endl;
		output << "c=1" << endl << endl;
		
		output << "######## calculate each molecules. need modify input file ##########" << endl << endl;
		output << "nmol=2" << endl << endl;
		output << "echo $nmol > mo_all.txt" << endl;
		output << "echo $nmol > nat_all.txt" << endl;
		output << "echo $nmol > nb_all.txt" << endl << endl;
		
		output << "if [ -f molecule_" << mol_label[mol_n1] << ".out ]; then" << endl;
		output << "	cat mo_moner_" << mol_label[mol_n1] << ".txt >> mo_all.txt" << endl;
		output << "	cat nat_" << mol_label[mol_n1] << ".txt >> nat_all.txt" << endl;
		output << "	cat nb_" << mol_label[mol_n1] << ".txt >> nb_all.txt" << endl;
		output << "else" << endl;
		output << "	" << zindo_folder << "/molecule/zindo1 <molecule_" << mol_label[mol_n1] << ".inp >molecule_" << mol_label[mol_n1] << ".out" << endl;
		output << "	cp mo_moner.txt mo_moner_" << mol_label[mol_n1] << ".txt >> mo_all.txt" << endl;
		output << "	cp nat.txt nat_" << mol_label[mol_n1] << ".txt >> nat_all.txt" << endl;
		output << "	cp nb.txt nb_" << mol_label[mol_n1] << ".txt >> nb_all.txt" << endl;		
		output << "	cat mo_moner_" << mol_label[mol_n1] << ".txt >> mo_all.txt" << endl;
		output << "	cat nat_" << mol_label[mol_n1] << ".txt >> nat_all.txt" << endl;
		output << "	cat nb_" << mol_label[mol_n1] << ".txt >> nb_all.txt" << endl;
		output << "fi" << endl << endl;

		output << "if [ -f molecule_" << mol_label[mol_n2] << ".out ]; then" << endl;
		output << "	cat mo_moner_" << mol_label[mol_n2] << ".txt >> mo_all.txt" << endl;
		output << "	cat nat_" << mol_label[mol_n2] << ".txt >> nat_all.txt" << endl;
		output << "	cat nb_" << mol_label[mol_n2] << ".txt >> nb_all.txt" << endl;
		output << "else" << endl;
		output << "	" << zindo_folder << "/molecule/zindo1 <molecule_" << mol_label[mol_n2] << ".inp >molecule_" << mol_label[mol_n2] << ".out" << endl;
		output << "	cp mo_moner.txt mo_moner_" << mol_label[mol_n2] << ".txt >> mo_all.txt" << endl;
		output << "	cp nat.txt nat_" << mol_label[mol_n2] << ".txt >> nat_all.txt" << endl;
		output << "	cp nb.txt nb_" << mol_label[mol_n2] << ".txt >> nb_all.txt" << endl;		
		output << "	cat mo_moner_" << mol_label[mol_n2] << ".txt >> mo_all.txt" << endl;
		output << "	cat nat_" << mol_label[mol_n2] << ".txt >> nat_all.txt" << endl;
		output << "	cat nb_" << mol_label[mol_n2] << ".txt >> nb_all.txt" << endl;
		output << "fi" << endl << endl;
		
		if (sign){
			output << "H_1=`sed -n " << coeff_H_lign << "p molecule_" << mol_label[mol_n1] << ".out | awk '{ print $" << coeff_H_row << " }'`" << endl;
			output << "L_1=`sed -n " << coeff_L_lign << "p molecule_" << mol_label[mol_n1] << ".out | awk '{ print $" << coeff_L_row << " }'`" << endl;
			output << "H_2=`sed -n " << coeff_H_lign << "p molecule_" << mol_label[mol_n2] << ".out | awk '{ print $" << coeff_H_row << " }'`" << endl;
			output << "L_2=`sed -n " << coeff_L_lign << "p molecule_" << mol_label[mol_n2] << ".out | awk '{ print $" << coeff_L_row << " }'`" << endl;
		}
		else{
			output << "H_1='1.0'" << endl;
			output << "L_1='1.0'" << endl;
			output << "H_2='1.0'" << endl;
			output << "L_2='1.0'" << endl << endl;
		}
		
		output << "####### finsih calculation for each molecule ###########" << endl << endl;
		output << "count=0" << endl << endl;
		output << "let \"mol_total=$nmol\"" << endl << endl;
		output << "let \"nelement_total=$mol_total*($mol_total-1)/2\"" << endl << endl;
		output << "echo 'NO of molecule:' $mol_tot" << endl;
		output << "echo $mol_total $nelement_total > split.out" << endl << endl;
		output << "for ((i=1;i<=mol_total;i++))" << endl;
		output << "do" << endl;
		output << "	for ((j=1;j<=i-1;j++))" << endl;
		output << "	do" << endl;
		output << "		let \"count=$count+1\"" << endl << endl;
		output << "		echo $count $nelement_total" << endl;
		output << "		echo $i $j >> split.out" << endl;
		output << "		echo $i $j > makeinp.inp" << endl;
		output << "		" << zindo_folder << "/split/makeinp" << endl << endl;
		output << "		echo ' ' >> temp2" << endl;
		output << "		echo ' $DATAIN' >> temp2" << endl;
		output << "		echo ' ' >> temp2" << endl << endl;
		output << "		echo $nlaya $nlayb $nlayc $a $b $c $nmol $i > translate.inp" << endl;
		output << "		" << zindo_folder << "/split/translate >> temp2" << endl;
		output << "		mv mo.tmp mo.txt" << endl << endl;
		output << "		echo $nlaya $nlayb $nlayc $a $b $c $nmol $j > translate.inp" << endl;
		output << "		" << zindo_folder << "/split/translate >> temp2" << endl;
		output << "		cat mo.tmp >> mo.txt" << endl << endl;
		output << "		echo ' ' >> temp2" << endl;
		output << "		echo ' $END' >> temp2" << endl << endl;
		output << "		" << zindo_folder << "/split/zindo-split-mod" << endl;
		output << "		cat split.tmp >> split.out" << endl << endl;
		output << "	done" << endl;
		output << "done" << endl << endl;
		output << "if [[ `wc -l split.out | awk '{print $1}'` -le 4 ]]; then" << endl;
		output << "	echo 'Error with file " << output_folder << "/frame_" << frame << "/dimer_"  <<  mol_label[mol_n1] << "_" << mol_label[mol_n2] << ".out" << "' >> " << log_file << ".log" << endl;
		output << "	J_H='0.0'" << endl;
		output << "	J_L='0.0'" << endl;
		output << "else" << endl;
		output << "	J_H=`grep ' i(1)" << setw(13) << n_electrons[mol_n1]/2 << " j(2)" << setw(13) << n_electrons[mol_n2]/2 << "' split.out | awk '{print $6}'`" << endl;
		output << "	J_L=`grep ' i(1)" << setw(13) << (n_electrons[mol_n1]/2)+1 << " j(2)" << setw(13) << (n_electrons[mol_n2]/2)+1 << "' split.out | awk '{print $6}'`" << endl;
		output << "fi" << endl << endl;
		output << "mv split.out dimer_"  <<  mol_label[mol_n1] << "_" << mol_label[mol_n2] << ".out" << endl;
		output << "echo " << frame << " " << mol_label[mol_n1] << " " << mol_label[mol_n2] << " $J_H $H_1 $H_2 $J_L $L_1 $L_2 >> frame_" << frame << ".out" << endl;
		output.close();
	}
	else 
		cerr << "Error opening " << output_filename.str().c_str() << endl;	
}

void Write_ZINDO_Files_MT(string input_file, string output_folder, string log_file, string zindo_folder){
	
	#pragma omp parallel for
	
	// Start calculation
	for (int i=0; i<n_frame; i++){
		double **mol1_cart, **mol2_cart;
		double **mol2_frac;
		
		Write_morange(input_file, i);
		
		for (int ii=0; ii<n_mol; ii++){
			
			mol1_cart = new double*[n_atom[ii]];
			
			for (int iii=0; iii<n_atom[ii]; iii++){	
				mol1_cart[iii] = new double[3];
				mol1_cart[iii][0] = x_cart[i][ii][iii];
				mol1_cart[iii][1] = y_cart[i][ii][iii];
				mol1_cart[iii][2] = z_cart[i][ii][iii];
			}	
					
			Write_INP(input_file, mol1_cart, i, ii);
			
			for (int iii=0; iii<n_atom[ii]; iii++){	
				delete [] mol1_cart[iii];
			}
			delete [] mol1_cart;
			
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
				
				int kk;
				for (int xx=0; xx<n_mol; xx++){
					if (mol_label[xx] == neigh_label[i][ii][jj]){
						kk = xx;
					}
				}
				
				mol1_cart = new double*[n_atom[ii]];
				mol2_cart = new double*[n_atom[kk]];
				mol2_frac = new double*[n_atom[kk]];
				
				for (int iii=0; iii<n_atom[ii]; iii++){	
					mol1_cart[iii] = new double[3];
					mol1_cart[iii][0] = x_cart[i][ii][iii];
					mol1_cart[iii][1] = y_cart[i][ii][iii];
					mol1_cart[iii][2] = z_cart[i][ii][iii];
					
				}
				for (int jjj=0; jjj<n_atom[jj+(ii+1)]; jjj++){	
					mol2_cart[jjj] = new double[3];
					mol2_cart[jjj][0] = x_cart[i][kk][jjj];
					mol2_cart[jjj][1] = y_cart[i][kk][jjj];
					mol2_cart[jjj][2] = z_cart[i][kk][jjj];
					
					mol2_frac[jjj] = new double[3];
				}
				
				//if(i==0){
				//	cout << i << " " << ii << " " << jj << endl;
				//	cout << mol2_cart[0][0] << " " << mol2_cart[0][1] << " " << mol2_cart[0][2] << endl;
				//}
				
				Mol_Cart_To_Frac(mol2_cart, mol2_frac, i, n_atom[kk]);
				for (int jjj=0; jjj<n_atom[kk]; jjj++){
					mol2_frac[jjj][0] = mol2_frac[jjj][0] + neigh_jump_vec_a[i][ii][jj];
					mol2_frac[jjj][1] = mol2_frac[jjj][1] + neigh_jump_vec_b[i][ii][jj];
					mol2_frac[jjj][2] = mol2_frac[jjj][2] + neigh_jump_vec_c[i][ii][jj];
				}
				Mol_Frac_To_Cart(mol2_frac, mol2_cart, i, n_atom[kk]);
				
				//if(i==0)
				//	cout << mol2_cart[0][0] << " " << mol2_cart[0][1] << " " << mol2_cart[0][2] << endl;
				
				Write_DAT(input_file, mol1_cart, mol2_cart, i, ii, kk);
				Write_CMD(input_file, zindo_folder, output_folder, log_file, i, ii, kk);
				
				for (int iii=0; iii<n_atom[ii]; iii++){	
					delete [] mol1_cart[iii];
				}
				
				for (int jjj=0; jjj<n_atom[kk]; jjj++){	
					delete [] mol2_cart[jjj];
					delete [] mol2_frac[jjj];
				}
				delete [] mol1_cart; delete [] mol2_cart; delete [] mol2_frac;
			}
		}
	}
	
	#pragma omp barrier
}

void Write_NB(string output_file, string input_folder){
	string tmp;
	stringstream file_nb;
	file_nb << input_folder << "/" << output_file << ".nb";
	
	FILE * pFile;

	pFile = fopen (file_nb.str().c_str(), "w");
	
	fprintf(pFile, "%d %d\n", n_frame, n_mol);
	for (int i=0; i<n_frame; i++){
		fprintf(pFile, "frame %d\n", i);
		for (int ii=0; ii<n_mol; ii++){
			fprintf(pFile, "molecule %d %d\n", mol_label[ii], int(neigh_label[i][ii].size()));
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
				fprintf(pFile, "%d\n", neigh_label[i][ii][jj]);
			}
		}	
	}
	
	fclose (pFile);
	
}

void Write_MC_YO(string input_file, string result_folder){
	string tmp;
	stringstream file_mc;
	file_mc << result_folder << "/" << input_file << ".mc_yo";
	
	FILE * pFile;

	pFile = fopen (file_mc.str().c_str(), "w");
	
	fprintf(pFile, "%d\n", n_frame);
	fprintf(pFile, "dt\n");
	fprintf(pFile, "0 0\n");
	fprintf(pFile, "lambda_i lambda_s T h_omega dist n_charges\n");
	fprintf(pFile, "F_norm\n");
	fprintf(pFile, "F_x F_y F_z\n");
	fprintf(pFile, "%d\n\n", n_mol);
	for (int i=0; i<n_frame; i++){
		fprintf(pFile, "N_frame=%d\n", i);
		fprintf(pFile, "%.4f %.4f %.4f %.4f %.4f %.4f\n", a[i], b[i], c[i], alpha_deg[i], beta_deg[i], gamma_deg[i]);
		for (int ii=0; ii<n_mol; ii++){
			fprintf(pFile, "0 %d 0\n", int(neigh_label[i][ii].size()));
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
				fprintf(pFile, "%d %e %e\n", neigh_label[i][ii][jj]+1, J_H[i][ii][jj]/1000.0, J_L[i][ii][jj]/1000.0);
			}
			fprintf(pFile, "%e %e %e 0\n", CM_x[i][ii], CM_y[i][ii], CM_z[i][ii]);
		}	
	}
	
	fclose (pFile);
	
}

void Write_MC(string input_file, string result_folder){
	string tmp;
	stringstream file_mc;
	file_mc << result_folder << "/" << input_file << ".mc";
	
	FILE * pFile;

	pFile = fopen (file_mc.str().c_str(), "w");
	
	fprintf(pFile, "%d %d\n", n_frame, n_mol);
	fprintf(pFile, "dt\n");
	fprintf(pFile, "lambda_i lambda_s T h_omega dist n_charges\n");
	fprintf(pFile, "F_norm\n\n");
	for (int i=0; i<n_frame; i++){
		fprintf(pFile, "frame %d\n", i);
		for (int ii=0; ii<n_mol; ii++){
			fprintf(pFile, "mol %d %d\n", mol_label[ii], int(neigh_label[i][ii].size()));
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
				fprintf(pFile, "%d %e %e\n", neigh_label[i][ii][jj], J_H[i][ii][jj]/1000.0, J_L[i][ii][jj]/1000.0);
			}
		}	
	}
	
	fclose (pFile);
	
}

int main(int argc, char **argv){

	int s;
	string input_file;
	string input_folder;
	string output_folder;
	string log_file;
	string zindo_folder;
	string result_folder;
	
	bool zindo = false;
	bool mc = false;
	
  	while ((s = getopt_long (argc, argv, "I:i:o:L:z:r:t:", NULL, NULL)) != -1){
      	switch (s){
			case 'I':
				input_file = optarg;
				break;
			case 'i':
				input_folder = optarg;
	  			break;
			case 'o':
				output_folder = optarg;
	  			break;
			case 'l':
				log_file = optarg;
	  			break;
			case 'z':
				zindo_folder = optarg;
	  			break;
			case 'r':
				result_folder = optarg;
	  			break;
	  		case 't':
	  			string a;
	  			a = optarg;
	  			if (a.compare("zindo") == 0)
	  				zindo = true;
	  			else if (a.compare("mc") == 0)
	  				mc = true;
		}
	}
	
	if (zindo){
		Read_XYZ(input_file, input_folder, false);
	}
	
	//if (mc){
	//	Read_NB(input_file, input_folder, false);
	//}
	
	Read_CELL(input_file, input_folder, false);
	Read_CM(input_file, input_folder, false);
	Read_ZIN(input_file, input_folder, false);
	
	if (zindo){
		Find_Neighbors_Sphere_MT(input_file, output_folder, false);
		Write_ZINDO_Files_MT(input_file, output_folder, log_file, zindo_folder);
		Write_NB(input_file, input_folder);
	}
	
	if (mc){
		Read_FULL(input_file, result_folder);
		Write_MC_YO(input_file, result_folder);
		Write_MC(input_file, result_folder);
	}

//==============================================================================

	if (zindo){
		for(int i=0; i<n_frame; i++){
			for(int ii=0; ii<n_mol; ii++){
				delete [] x_cart[i][ii];
				delete [] y_cart[i][ii];
				delete [] z_cart[i][ii];
				delete [] atomic_number[i][ii];
				delete [] atomic_mass[i][ii];
				delete [] atomic_valence[i][ii];
				delete [] symbol[i][ii];
			}
			
			delete [] x_cart[i];
			delete [] y_cart[i];
			delete [] z_cart[i];
			delete [] atomic_number[i];
			delete [] atomic_mass[i];
			delete [] atomic_valence[i];
			delete [] symbol[i];
		}
		delete [] x_cart; delete [] y_cart; delete [] z_cart; 
		delete [] atomic_number; delete [] atomic_mass; delete [] atomic_valence; delete [] symbol; 
	}
	
	if (mc){
		J_H.clear(); J_L.clear(); 
	}
	
	for(int i=0; i<n_frame; i++){
		delete [] CM_x[i];
		delete [] CM_y[i];
		delete [] CM_z[i];
	}
	delete [] mol_label; delete [] J;
	delete [] CM_x; delete [] CM_y; delete [] CM_z; delete [] n_electrons;
	
	delete [] n_atom;
	
	delete [] a; delete [] b; delete [] c; delete [] alpha_deg; delete [] beta_deg; delete [] gamma_deg;
	delete [] temp_alpha_cos; delete [] temp_beta_sin; delete [] temp_beta_cos; delete [] temp_gamma_sin;
	delete [] temp_gamma_cos; delete [] temp_beta_term; delete [] temp_gamma_term;
	
	neigh_label.clear();
	neigh_jump_vec_a.clear(); neigh_jump_vec_b.clear(); neigh_jump_vec_c.clear();
	
	return 0;
}

