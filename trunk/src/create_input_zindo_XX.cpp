#include <iostream> //Entrées-sorties standard
#include <fstream> //Entrées-sorties sur fichiers
#include <string> //Chaines de caracteres
//#include <limits> //Pour aller à la fin d'une ligne, par exemple
#include <iomanip> //Manipulation des sorties, par exemple format scientifique
#include <sstream> //Pour les conversion entre types de variables
#include <math.h>
#include <getopt.h>

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

int ****displ_vec; bool ***neighbors;

void Read_XYZ(string input_file, bool print_results){
	string tmp;
	stringstream file_xyz;
	file_xyz << input_file << ".xyz";
	
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
		mol_label = new int[n_mol];
		symbol = new char**[n_frame];
		
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

void Read_CELL(string input_file, bool print_results){
	string tmp;
	stringstream file_cell;
	file_cell << input_file << ".cell";
	
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

void Read_CM(string input_file, bool print_results){
	string tmp;
	stringstream file_cm;
	file_cm << input_file << ".cm";
	
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

double* Cartesian_To_Fractional(double* Dist_Cart, int i){
	double *tmp;
	tmp = new double[3];
	
	double x = Dist_Cart[0];
	double y = Dist_Cart[1];
	double z = Dist_Cart[2];
	
	z = (z/temp_gamma_term[i]) / c[i];
	y = ((y-z*c[i]*temp_beta_term[i])/temp_gamma_sin[i]) / b[i];
	x = (x-y*b[i]*temp_gamma_cos[i]-z*c[i]*temp_beta_cos[i]) / a[i];
	tmp[0] = x; tmp[1] = y; tmp[2] = z;
	
	return tmp;
	delete [] tmp;
}

double* Fractional_To_Cartesian(double* Dist_Frac, int i){
	double *tmp;
	tmp = new double[3];
	
	double x = Dist_Frac[0];
	double y = Dist_Frac[1];
	double z = Dist_Frac[2];
		
	x = x*a[i] + y*b[i]*temp_gamma_cos[i] + z*c[i]*temp_beta_cos[i];
	y = y*b[i]*temp_gamma_sin[i] + z*c[i]*temp_beta_term[i];
	z = z*c[i]*temp_gamma_term[i];
	tmp[0] = x; tmp[1] = y; tmp[2] = z;
	
	return tmp;
	delete [] tmp;
}

double** Mol_Frac_To_Cart(double** Coord_Frac, int frame, int atom){
	double **Coord_Cart;
	
	Coord_Cart = new double*[atom];
	for (int i=0; i<atom; i++){
		Coord_Cart[i] = new double[3];
		Coord_Cart[i] = Fractional_To_Cartesian(Coord_Frac[i], frame);
	}
	
	return Coord_Cart;

	for (int i=0; i<atom; i++){
		delete [] Coord_Cart[i];
	}
	delete [] Coord_Cart;
}

double** Mol_Cart_To_Frac(double** Coord_Cart, int frame, int atom){
	double **Coord_Frac;
	
	Coord_Frac = new double*[atom];
	for (int i=0; i<atom; i++){
		Coord_Frac[i] = new double[3];
		Coord_Frac[i] = Cartesian_To_Fractional(Coord_Cart[i], frame);
	}
	
	return Coord_Frac;

	for (int i=0; i<atom; i++){
		delete [] Coord_Frac[i];
	}
	delete [] Coord_Frac;
}

void Find_Neighbors_Sphere(bool print_results){

	displ_vec = new int***[n_frame];
	neighbors = new bool**[n_frame];
	
	for(int i=0; i<n_frame; i++){
		displ_vec[i] = new int**[n_mol];
		neighbors[i] = new bool*[n_mol];
		
		for(int ii=0; ii<n_mol; ii++){
			displ_vec[i][ii] = new int*[n_mol];
			neighbors[i][ii] = new bool[n_mol];
			
			for(int jj=0; jj<n_mol; jj++){
				displ_vec[i][ii][jj] = new int[3];
				neighbors[i][ii][jj] = 0;
				
				for(int k=0; k<3; k++){
					displ_vec[i][ii][jj][k] = 0;
				}
			}
		}
	}

	double CutOff_square = pow(cutoff, 2);
	
	double Dist_Norm_square = 0.0;
	double *Dist_Cart, *Dist_Frac;
	Dist_Cart = new double[3];
	Dist_Frac = new double[3];
	
	for (int i=0; i<n_frame; i++){
		for (int ii=0; ii<n_mol; ii++){
			if(J[ii]){
				for (int jj=ii+1; jj<n_mol; jj++){
					Dist_Cart[0] = CM_x[i][ii] - CM_x[i][jj];
					Dist_Cart[1] = CM_y[i][ii] - CM_y[i][jj];
					Dist_Cart[2] = CM_z[i][ii] - CM_z[i][jj];
					
					Dist_Frac = Cartesian_To_Fractional(Dist_Cart, i);

					// Periodic boundary conditions calculations
					for (int k=0; k<3; k++){
						if (fabs(Dist_Frac[k]) > 0.5 && pbc[k]){
							if (Dist_Frac[k] < 0.0){
								Dist_Frac[k] = 1.0 + Dist_Frac[k];
								displ_vec[i][ii][jj][k] = 1;
							}
							else{
								Dist_Frac[k] = 1.0 - Dist_Frac[k];
								displ_vec[i][ii][jj][k] = -1;
							}
						}
					}
					
					Dist_Cart = Fractional_To_Cartesian(Dist_Frac, i);
					
					Dist_Norm_square = pow(Dist_Cart[0],2) + pow(Dist_Cart[1],2) + pow(Dist_Cart[2],2);
					
					if (Dist_Norm_square < CutOff_square && Dist_Norm_square != 0){
						neighbors[i][ii][jj] = 1;
					}
				}
			}
		}
	}
	
	int tmp;
	if(print_results){
		for (int i=0; i<n_frame; i++){
			cout << "frame " << i << endl;
			for (int ii=0; ii<n_mol; ii++){
				tmp = 0;
				for (int jj=0; jj<n_mol; jj++){
					if (neighbors[i][ii][jj]){
						tmp++;
					}
				}
				cout << "Molecule " << mol_label[ii] << ": " << tmp << " Neighbors\n";
				for (int jj=0; jj<n_mol; jj++){
					if (neighbors[i][ii][jj]){
						cout << "Vector " << displ_vec[i][ii][jj][0] << " " << displ_vec[i][ii][jj][1] << " " << displ_vec[i][ii][jj][2] << endl;
					}
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
	ofstream output(output_filename.str().c_str(), ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier
	
	if (output){
		output << n_atom[mol_n1] << " " << n_atom[mol_n2] << endl;
		for (int iii=0; iii<n_atom[mol_n1]; iii++)
			output << fixed << setw(10) << setprecision(5) << mol1[iii][0] << setw(10) << mol1[iii][1] << setw(10) << mol1[iii][2] << setw(4) << atomic_number[frame][mol_n1][iii] << endl;
		
		for (int jjj=0; jjj<n_atom[mol_n2]; jjj++)
			output << fixed << setw(10) << setprecision(5) << mol2[jjj][0] << setw(10) << mol2[jjj][1] << setw(10) << mol2[jjj][2] << setw(4) << atomic_number[frame][mol_n2][jjj] << endl;
		
		output.close();
	}
	
	else 
		cerr << "Error opening " << output_filename.str().c_str() << endl;
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
	ofstream output(output_filename.str().c_str(), ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier
	
	if (output){
		output << " $TITLEI" << endl;
		output << "INPUT CREATED BY Mon cul c'est du poulet" << endl << endl;
		output << " $END" << endl << endl;
		output << " $CONTRL" << endl << endl;
		output << " SCFTYP     RHF      RUNTYP     ENERGY     ENTTYP     COORD    UNITS      ANGS" << endl << endl;
		output << "INTTYP          1   IAPX            3   NAT          " << setw(6) << n_atom[mol_n] << "   NEL       " << setw(4) <<  n_electrons[mol_n] << " " << endl;
		output << "MULT            1   ITMAX          50 " << endl << endl;
		output << " SCFTOL   0.0000010  " << endl << endl;
		output << "! ***** Basis set and C.I. size inFormation *****" << endl << endl;
		output << "   DYNAL(1) =    0" << setw(5) << n_s << setw(5) << n_p << setw(5) << n_d << "  0  4000" << setw(5) << n_electrons[mol_n]/2 << " " << endl << endl;
		output << "! ***** Interaction factors *****" << endl << endl;
		output << "   INTFA(1) =   1.00000  1.26700  0.58500  1.00000  1.00000  1.00000" << endl << endl;
		output << "! ***** Output file name *****" << endl << endl;
		output << "   ONAME = file99" << endl << endl;
		output << " $END" << endl << endl;
		output << " $DATAIN" << endl;
		for (int iii=0; iii<n_atom[mol_n]; iii++){
				output << fixed << setw(10) << setprecision(5) << mol[iii][0] << setw(10) << mol[iii][1] << setw(10) << mol[iii][2] << setw(4) << atomic_number[frame][mol_n][iii] << endl;
		}
		output << endl;
		output << " $END" << endl << endl;
		output << " $CIINPU" << endl << endl;
		output << "! ***** Configuration Interaction specification *****" << endl << endl;
		output << "    2    1 " << setw(2) << n_exc << "    1    0    0    0    1    1    1 " << setw(2) << n_exc << endl; // A faire
		output << "  -80000.0  0.000000" << endl << endl;
		output << "00001" << setw(5) << n_electrons[mol_n]/2  << setw(5) << n_electrons[mol_n]/2  << endl;
		output << "00021" << setw(5) << (n_electrons[mol_n]/2)+1  << setw(5) << n_electrons[mol_n]/2 << setw(5) << (n_electrons[mol_n]/2)+1 << setw(5) << n_electrons[mol_n]/2 << endl << endl << endl << endl;
		output << " $END";
	}
	else 
		cerr << "Error opening " << output_filename.str().c_str() << endl;
}

void Write_morange(string input_file, int frame){
	
	stringstream output_filename;
	stringstream s_frame;
	
	s_frame << frame;
	output_filename << input_file.c_str() << "/frame_" << s_frame.str().c_str() << "/morange.inp";
	ofstream output(output_filename.str().c_str(), ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier
	
	if (output){
		output << n_electrons[0]/2 - 4 << " " << n_electrons[0]/2 + 4 << endl;
		output << n_electrons[0]/2 - 4 << " " << n_electrons[0]/2 + 4;
		output.close();
	}
	else 
		cerr << "Error opening " << output_filename.str().c_str() << endl;
}

void Write_CMD(string input_file, string output_folder, int frame, int mol_n1, int mol_n2){

	stringstream output_filename;
	stringstream s_frame, s_mol_n1, s_mol_n2;
	
	s_frame << frame;
	s_mol_n1 << mol_label[mol_n1];
	s_mol_n2 << mol_label[mol_n2];
	output_filename << input_file.c_str() << "/frame_" << s_frame.str().c_str() << "/run_" << s_mol_n1.str().c_str() << "_" << s_mol_n2.str().c_str() << ".cmd";
	ofstream output(output_filename.str().c_str(), ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier
	
	if (output){
		output << "#!/bin/bash" << endl << endl;
		output << "cp cell_"  <<  mol_n1 << "_" << mol_n2 << ".dat cell.dat" << endl << endl;
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
		output << "/home/output/David/Aijun/zindo-split/molecule/zindo1 <molecule_" << mol_n1 << ".inp >molecule_" << mol_n1 << ".out" << endl;
		output << "cat mo_moner.txt >> mo_all.txt" << endl;
		output << "cat nat.txt >> nat_all.txt" << endl;
		output << "cat nb.txt >> nb_all.txt" << endl << endl;
		output << "/home/output/David/Aijun/zindo-split/molecule/zindo1 <molecule_" << mol_n2 << ".inp >molecule_" << mol_n2 << ".out" << endl;
		output << "cat mo_moner.txt >> mo_all.txt" << endl;
		output << "cat nat.txt >> nat_all.txt" << endl;
		output << "cat nb.txt >> nb_all.txt" << endl << endl;
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
		output << "		/home/output/David/Aijun/zindo-split/split/makeinp" << endl << endl;
		output << "		echo ' ' >> temp2" << endl;
		output << "		echo ' $DATAIN' >> temp2" << endl;
		output << "		echo ' ' >> temp2" << endl << endl;
		output << "		echo $nlaya $nlayb $nlayc $a $b $c $nmol $i > translate.inp" << endl;
		output << "		/home/output/David/Aijun/zindo-split/split/translate >> temp2" << endl;
		output << "		mv mo.tmp mo.txt" << endl << endl;
		output << "		echo $nlaya $nlayb $nlayc $a $b $c $nmol $j > translate.inp" << endl;
		output << "		/home/output/David/Aijun/zindo-split/split/translate >> temp2" << endl;
		output << "		cat mo.tmp >> mo.txt" << endl << endl;
		output << "		echo ' ' >> temp2" << endl;
		output << "		echo ' $END' >> temp2" << endl << endl;
		output << "		/home/output/David/Aijun/zindo-split/split/zindo-split-mod2" << endl;
		output << "		cat split.tmp >> split.out" << endl << endl;
		output << "	done" << endl;
		output << "done" << endl;
		output << "mv split.out " << output_folder << "/frame_" << frame << "/dimer_"  <<  mol_n1 << "_" << mol_n2 << ".out" << endl;
		output << "rm *.out" << endl;
	}
	else 
		cerr << "Error opening " << output_filename.str().c_str() << endl;	
}

void Write_ZINDO_Files(string input_file, string output_folder){
	double **mol1_cart, **mol2_cart;
	double **mol2_frac;
	
	for (int i=0; i<n_frame; i++){
		
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
			
			for (int jj=ii+1; jj<n_mol; jj++){
				if (neighbors[i][ii][jj]){
					mol1_cart = new double*[n_atom[ii]];
					mol2_cart = new double*[n_atom[jj]];
					mol2_frac = new double*[n_atom[jj]];
					
					for (int iii=0; iii<n_atom[ii]; iii++){	
						mol1_cart[iii] = new double[3];
						mol1_cart[iii][0] = x_cart[i][ii][iii];
						mol1_cart[iii][1] = y_cart[i][ii][iii];
						mol1_cart[iii][2] = z_cart[i][ii][iii];
						
					}
					for (int jjj=0; jjj<n_atom[jj]; jjj++){	
						mol2_cart[jjj] = new double[3];
						mol2_cart[jjj][0] = x_cart[i][jj][jjj];
						mol2_cart[jjj][1] = y_cart[i][jj][jjj];
						mol2_cart[jjj][2] = z_cart[i][jj][jjj];
						
						mol2_frac[jjj] = new double[3];
					}
					
					//if(i==0){
					//	cout << i << " " << ii << " " << jj << endl;
					//	cout << mol2_cart[0][0] << " " << mol2_cart[0][1] << " " << mol2_cart[0][2] << endl;
					//}
					
					mol2_frac = Mol_Cart_To_Frac(mol2_cart, i, n_atom[jj]);
					for (int jjj=0; jjj<n_atom[jj]; jjj++){
						for (int k=0; k<3; k++){
							mol2_frac[jjj][k] = mol2_frac[jjj][k] - displ_vec[i][ii][jj][k];
						}
					}
					mol2_cart = Mol_Frac_To_Cart(mol2_frac, i, n_atom[jj]);
					
					//if(i==0)
					//	cout << mol2_cart[0][0] << " " << mol2_cart[0][1] << " " << mol2_cart[0][2] << endl;
					
					Write_DAT(input_file, mol1_cart, mol2_cart, i, ii, jj);
					Write_CMD(input_file, output_folder, i, ii, jj);
					
					for (int iii=0; iii<n_atom[ii]; iii++){	
						delete [] mol1_cart[iii];
					}
					
					for (int jjj=0; jjj<n_atom[jj]; jjj++){	
						delete [] mol2_cart[jjj];
						delete [] mol2_frac[jjj];
					}
					delete [] mol1_cart; delete [] mol2_cart; delete [] mol2_frac; 		
					
				}
			}
		}
	}		
}

int main(int argc, char **argv){

	int s;
	string input_file;
	string output_folder;
	
  	while ((s = getopt_long (argc, argv, "i:o:", NULL, NULL)) != -1){
      	switch (s){
			case 'i':
				input_file = optarg;
	  			break;
			case 'o':
				output_folder = optarg;
	  			break;
		}
	}
	
	Read_XYZ(input_file, false);
	Read_CELL(input_file, false);
	Read_CM(input_file, false);
	Find_Neighbors_Sphere(false);
	Write_ZINDO_Files(input_file, output_folder);

	for(int i=0; i<n_frame; i++){
		for(int ii=0; ii<n_mol; ii++){
			for(int jj=0; jj<n_mol; jj++){
				delete [] displ_vec[i][ii][jj];
			}
			delete [] displ_vec[i][ii];
			delete [] neighbors[i][ii];
		}
		delete [] displ_vec[i];
		delete [] neighbors[i];
	}
	delete [] displ_vec; 
	delete [] neighbors;

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
		delete [] CM_x[i];
		delete [] CM_y[i];
		delete [] CM_z[i];
	}
	delete [] x_cart; delete [] y_cart; delete [] z_cart; 
	delete [] atomic_number; delete [] atomic_mass; delete [] atomic_valence;
	delete [] symbol; delete [] mol_label; delete [] J;
	delete [] CM_x; delete [] CM_y; delete [] CM_z; delete [] n_electrons;
	
	delete [] n_atom;
	
	delete [] a; delete [] b; delete [] c; delete [] alpha_deg; delete [] beta_deg; delete [] gamma_deg;
	delete [] temp_alpha_cos; delete [] temp_beta_sin; delete [] temp_beta_cos; delete [] temp_gamma_sin;
	delete [] temp_gamma_cos; delete [] temp_beta_term; delete [] temp_gamma_term;
	
	return 0;
}
