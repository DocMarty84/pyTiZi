#include <iostream> //Entrées-sorties standard
#include <fstream> //Entrées-sorties sur fichiers
#include <string> //Chaines de caracteres
//#include <limits> //Pour aller à la fin d'une ligne, par exemple
#include <iomanip> //Manipulation des sorties, par exemple format scientifique
#include <sstream> //Pour les conversion entre types de variables
#include <math.h>
#include <getopt.h>
#include <stdio.h>

using namespace std;

double **CM_x, **CM_y, **CM_z; // Center of masses
int n_frame, n_mol, *mol_label; // Number of frame, molecule and label of molecules

double *a, *b, *c, *alpha_deg, *beta_deg, *gamma_deg; // Cell parameters

bool ***neighbors; int **n_neighbors; int **mol_label_initial, **mol_label_final;

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

		input >> tmp >> tmp >> tmp >> tmp;
		input >> tmp >> tmp;
		for(int i=0; i<n_frame; i++){
			input >> tmp >> tmp;
			input >> a[i] >> b[i] >> c[i] >> alpha_deg[i] >> beta_deg[i] >> gamma_deg[i];
			input >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp;
		}
		
		input.close();
		
		if (print_results){
			for(int i=0; i<n_frame; i++){
				cout << "frame " << i << endl;
				cout << a[i] << " " << b[i] << " " << c[i] << " " << alpha_deg[i] << " " << beta_deg[i] << " " << gamma_deg[i] << endl;
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
		mol_label = new int[n_mol];
		
		for(int i=0; i<n_frame; i++){
			CM_x[i] = new double[n_mol];
			CM_y[i] = new double[n_mol];
			CM_z[i] = new double[n_mol];
		}
		
		for(int i=0; i<n_frame; i++){
			input >> tmp >> tmp;
			for(int ii=0; ii<n_mol; ii++){
				input >> tmp >> mol_label[ii] >> tmp >> CM_x[i][ii] >> CM_y[i][ii] >> CM_z[i][ii] >> tmp;
			}
		}
		
		input.close();
		
		if(print_results){
			for(int i=0; i<n_frame; i++){
				cout << "frame " << i << endl;
				for(int ii=0; ii<n_mol; ii++){
					cout << "molecule " << mol_label[ii] << " " << tmp << " " << CM_x[i][ii] << " " << CM_y[i][ii] << " " << CM_z[i][ii] << " " << tmp << endl;
				}
			}
		}
	}
	else{
		cerr << "Error opening " << file_cm.str().c_str() << endl;
		exit(1);
	}
}

void Read_NB(string input_file, bool print_results){
	string tmp;
	stringstream file_cell;
	file_cell << input_file << ".nb";
	
	ifstream input(file_cell.str().c_str(), ios::in);
	if (input){
		input >> n_frame >> n_mol;

		for (int i=0; i<n_frame; i++){
			input >> tmp >> tmp;
			for (int ii=0; ii<n_mol; ii++){
				input >> tmp >> 
			}
		}
	}
		
		input.close();
		
		if (print_results){
		}
		
	}
	else{
		cerr << "Error opening " << file_cell.str().c_str() << endl;
		exit(1);
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
	
	Read_CELL(input_file, false);
	Read_CM(input_file, false);

	for(int i=0; i<n_frame; i++){
		for(int ii=0; ii<n_mol; ii++){
			delete [] neighbors[i][ii];
		}
		delete [] neighbors[i];
		delete [] n_neighbors[i];
	}
	delete [] neighbors;
	delete [] n_neighbors;

	for(int i=0; i<n_frame; i++){
		delete [] CM_x[i];
		delete [] CM_y[i];
		delete [] CM_z[i];
	}
	
	delete [] a; delete [] b; delete [] c; delete [] alpha_deg; delete [] beta_deg; delete [] gamma_deg;
	
	return 0;
}
