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
#include "/home/nicolas/Downloads/alglib_cpp/out/matinv.h"

using namespace std;

int n_mo;
int homo_1=0, homo_2=0, lumo_1=0, lumo_2=0;
vector<double> E;
vector< vector<double> > H_KS, S, C, C_inv, I;

void Print_Coeff(vector< vector<double> > Coeff, string title){

	cout << title << endl << endl;

	int n_mo_rounded = (n_mo - (n_mo % 10))/10;
	
	for (int count_mo=0; count_mo<n_mo_rounded; count_mo++) {
		
		for(int j=0; j<11; j++) {

			if(j==1)
				printf(" MOs :%7d", j+(count_mo*10));
			
			else if(j>1)
				printf("%9d", j+(count_mo*10));
		}
		printf("\n");
		
		for(int i=0; i<n_mo; i++) {
			for(int j=0; j<11; j++) {

				if(j==1)
					printf("%6d%9.4lf", i+1, Coeff[i][(j-1)+(count_mo*10)]);
				
				else if(j>1)
					printf("%9.4lf", Coeff[i][(j-1)+(count_mo*10)]);
			}
			printf("\n");
		}
		printf("\n");
	}
	
	if (n_mo % 10 != 0) {
		
		for(int j=0; j<(n_mo % 10)+1; j++) {

			if(j==1)
				printf(" MOs :%7d", j+(n_mo_rounded*10));
			
			else if(j>1)
				printf("%9d", j+(n_mo_rounded*10));
		}
		printf("\n");
		
		for(int i=0; i<n_mo; i++) {
			for(int j=0; j<(n_mo % 10)+1; j++) {
				
				if(j==1)
					printf("%6d%9.4lf", i+1, Coeff[i][(j-1)+(n_mo_rounded*10)]);
				
				else if(j>1)
					printf("%9.4f", Coeff[i][(j-1)+(n_mo_rounded*10)]);
			}
			printf("\n");
		}
	}
	printf("\n");
	
	cout << "=====================================" << endl << endl;
}

void Print_Triangular(vector< vector<double> > Triang, string title){

	cout << title << endl << endl;

	int n_mo_rounded = (n_mo - (n_mo % 4))/4;
	
	for (int count_mo=0; count_mo<n_mo_rounded; count_mo++) {

		int k = 1;
		
		for(int j=0; j<5; j++) {

			if(j==1){
				//printf(" column%12d", j+(count_mo*4));
				cout << " column" << setw(12) << j+(count_mo*4);
			}
			
			else if(j>1){
				//printf("%22d", j+(count_mo*4));
				cout << setw(22) << j+(count_mo*4);
			}
		}
		//printf("\n row\n");
		cout << endl << " row" << endl;
		
		for(int i=count_mo*4; i<n_mo; i++) {
			
			if (k<5) {
			  k++;
			}

			for (int j=0; j<k; j++) {

				if(j==1){
					//printf("%5d%24.14g", i+1, S[i][(j-1)+(count_mo*4)]);
					cout << setw(5) << i+1 << setprecision(14) << setw(24) << scientific << Triang[i][(j-1)+(count_mo*4)];
				}
				
				else if(j>1){
					//printf("%22.14g", S[i][(j-1)+(count_mo*4)]);
					cout << setprecision(14) << setw(22) << scientific << Triang[i][(j-1)+(count_mo*4)];
				}	
			}
			//printf("\n");
			cout << endl;
		}
		//printf("\n");
		cout << endl;
	}
	
	if (n_mo % 4 != 0) {

		int k = 1;
		
		for(int j=0; j<5; j++) {

			if(j==1){
				//printf(" column%12d", j+(n_mo_rounded*4));
				cout << " column" << setw(12) <<  j+(n_mo_rounded*4);
			}
			
			else if(j>1){
				//printf("%22d", j+(n_mo_rounded*4));
				cout << setw(22) << j+(n_mo_rounded*4);
			}
		}
		//printf("\n row\n");
		cout << endl << " row" << endl;
		
		for(int i=n_mo-(n_mo%4); i<n_mo; i++) {
			
			if (k<5) {
			  k++;
			}

			for (int j=0; j<k; j++) {

				if(j==1){
					//printf("%5d%24.14g", i+1, S[i][(j-1)+(n_mo_rounded*4)]);
					cout << setw(5) << i+1 << setprecision(14) << setw(24) << scientific << Triang[i][(j-1)+(n_mo_rounded*4)];
				}
				
				else if(j>1){
					//printf("%22.14g", S[i][(j-1)+(n_mo_rounded*4)]);
					cout << setprecision(14) << setw(22) << scientific << Triang[i][(j-1)+(n_mo_rounded*4)];
				}
			}
			cout << endl;
		}
	}
	cout << endl;
	
	cout << "=====================================" << endl << endl;
}

void Print_Table(vector<double> Tab, string title){
	
	cout << title << endl << endl;
	
	for (int i=0; i<n_mo; i++){
		cout << setw(5) << i+1 << setprecision(14) << setw(30) << scientific << Tab[i] << endl;
	}
	
	cout << endl;
	
	cout << "=====================================" << endl << endl;
}

vector< vector<double> > Matrix_Product (vector< vector<double> > M1, vector< vector<double> > M2) {

	vector< vector<double> > M_res;
	
	for(int i=0; i<n_mo; i++) {
		M_res.push_back( vector<double> ());
		
		for(int j=0; j<n_mo; j++){
			M_res[i].push_back(0.0);
		}
	}

	#pragma omp parallel for
	for (unsigned int i=0; i<M1[0].size(); i++) {
		double sum = 0.0;
		for (unsigned int j=0; j<M1[0].size(); j++) {
			sum = 0.0;
			for (unsigned int k=0; k<M1[0].size(); k++) {
				sum += M1[i][k] * M2[k][j];
			}
			//if (sum < 1e-13) {
			//	sum = 0.0;
			//}
			M_res[i][j] = sum;
		}
	}
	#pragma omp barrier
	
	return(M_res);

}

void Read_DIMER(string input_file, bool print_result){
	//stringstream file_dimer;
	//file_dimer << input_file << ".full";
	
	ifstream input(input_file.c_str(), ios::in);
	
	if (input){
		string str;
	
		while (!input.eof()){
			getline(input, str);
			//cout << str << endl;

			// Read number of MO
			if (str.find("Total nr. of (C)SFOs (summation over all irreps) :") != string::npos) {
				size_t length;
				char buffer[10];
				length = str.copy(buffer, 9, 51);
				buffer[length] = '\0';
				n_mo = atoi(buffer);
				
				if (n_mo%2 != 0){
					cerr << "[ERROR] Number of SFO odd (" << n_mo << ")! Unable to calculate J. Exiting..." << endl;
					S.clear(); C.clear(); C_inv.clear(); I.clear(); E.clear();
					exit(1); 
				}
			}
			
			// Read HOMO and LUMO numbers
			else if (str.find("SFO  (index         Fragment          Generating    Expansion in Fragment Orbitals") != string::npos) {

				do {
					getline(input, str);
				} 
				while (str.find("--------------------------------------------------------------------------------------") == string::npos);
				
				int count = 0;
				do {
					getline(input, str);
					if (str.find(" frag1 ") != string::npos) {
						count++;
					}
					if (str.find(" 2.000 ") != string::npos && str.find(" frag1 ") != string::npos) {
						homo_1++;
					}
					if (str.find(" 2.000 ") != string::npos && str.find(" frag2 ") != string::npos) {
						homo_2++;
					}
					getline(input, str);
				} 
				while (str.find("eV)") != string::npos);

				homo_1 = homo_1 - 1;
				homo_2 = homo_2 - 1 + count;
				
				lumo_1 = homo_1 + 1;
				lumo_2 = homo_2 + 1;
			    
			}
			
			// Read C matrix
			else if (str.find("MOs expanded in CFs+SFOs") != string::npos) {
				
				for(int i=0; i<n_mo; i++) {
					C.push_back( vector<double> ());
					C_inv.push_back( vector<double> ());
					
					for(int j=0; j<n_mo; j++){
						C[i].push_back(0.0);
						C_inv[i].push_back(0.0);
					}
				}
				
				// Specific variables for inversion library
				ap::real_2d_array a;
				a.setlength(n_mo, n_mo);
				
				int n_mo_rounded = (n_mo - (n_mo % 10))/10;
				double Z;
				
				for (int count_mo=0; count_mo<n_mo_rounded; count_mo++) {
					do {
						getline(input, str);
					} 
					while (str.find("CF+SFO") == string::npos);
					
					for(int i=0; i<n_mo; i++) {
						for(int j=0; j<11; j++) {
							//input >> str;
							input >> Z;

							if(j>0){
								//C[i][(j-1)+(count_mo*10)] = atof(str.c_str());
								//a(i,(j-1)+(count_mo*10)) = atof(str.c_str());
								C[i][(j-1)+(count_mo*10)] = Z;
								a(i,(j-1)+(count_mo*10)) = Z;
							}		
						}
					}
				}
				
				if (n_mo % 10 != 0) {
					do {
						getline(input, str);
					} 
					while (str.find("CF+SFO") == string::npos);
					
					for(int i=0; i<n_mo; i++) {
						for(int j=0; j<(n_mo % 10)+1; j++) {
							//input >> str;
							input >> Z;
	
							if(j>0){
								//C[i][(j-1)+(n_mo_rounded*10)] = atof(str.c_str());
								//a(i,(j-1)+(n_mo_rounded*10)) = atof(str.c_str());
								C[i][(j-1)+(n_mo_rounded*10)] = Z;
								a(i,(j-1)+(n_mo_rounded*10)) = Z;
							}		
						}
					}
				}


				// Calculates inverse matrix
				
				int info;
				matinvreport rep;
				rmatrixinverse(a, n_mo, info, rep);
				
				for(int i=0; i<n_mo; i++) {
					for(int j=0; j<n_mo; j++) {
						C_inv[i][j] = a(i,j);
					}
				}
			}

			// Read S matrix
			else if (str.find("SFO Overlap Matrix (valence part only)") != string::npos) {

				for(int i=0; i<n_mo; i++) {
					S.push_back( vector<double> ());
					
					for(int j=0; j<n_mo; j++){
						S[i].push_back(0.0);
					}
				}
				
				int n_mo_rounded = (n_mo - (n_mo % 4))/4;
				double Z;
				
				for (int count_mo=0; count_mo<n_mo_rounded; count_mo++) {
					do {
						getline(input, str);
					} 
					while (str.find("row") == string::npos);
					
					int k = 1;
					
					for(int i=count_mo*4; i<n_mo; i++) {
						
						if (k<5) {
						  k++;
						}

						for (int j=0; j<k; j++) {
							//input >> str;
							input >> Z;
							
							if(j>0){
								//S[i][(j-1)+(count_mo*4)] = atof(str.c_str());
								//S[(j-1)+(count_mo*4)][i] = atof(str.c_str());
								S[i][(j-1)+(count_mo*4)] = Z;
								S[(j-1)+(count_mo*4)][i] = Z;
							}	
						}
					}
				}
				
				if (n_mo % 4 != 0) {
					do {
						getline(input, str);
					} 
					while (str.find("row") == string::npos);
					
					int k = 1;
					
					for(int i=n_mo-(n_mo%4); i<n_mo; i++) {
						
						if (k<5) {
						  k++;
						}

						for (int j=0; j<k; j++) {
							//input >> str;
							input >> Z;
							
							if(j>0){
								//S[i][(j-1)+(n_mo_rounded*4)] = atof(str.c_str());
								//S[(j-1)+(n_mo_rounded*4)][i] = atof(str.c_str());
								S[i][(j-1)+(n_mo_rounded*4)] = Z;
								S[(j-1)+(n_mo_rounded*4)][i] = Z;
							}	
						}
					}
				}
			}
			
			// Read Energies
			else if (str.find("Orbital Energies, per Irrep and Spin:") != string::npos) {
			
				do {
					getline(input, str);
				} 
				while (str.find(" A") == string::npos);
				
				int mo;
				float occ;
				double E_au;
				
				for(int i=0; i<n_mo; i++) {
					input >> mo >> occ;
					
					if(occ == 2.0)
						input >> E_au >> str >> str;
					
					else
						input >> E_au >> str;
						
					E.push_back(E_au*27.211383);
					//E.push_back(E_au);
				}
			}
		}
		
		input.close();
		
	}
	
	else{
		cerr << "Error opening " << input_file.c_str() << endl;
		exit(1);			
	}
	
	if (print_result) {
		
		cout << "N_SFO = " << n_mo << endl;
		cout << "HOMO frag1 = " << homo_1+1 << endl;
		cout << "LUMO frag1 = " << lumo_1+1 << endl;
		cout << "HOMO frag2 = " << homo_2+1 << endl;
		cout << "LUMO frag2 = " << lumo_2+1 << endl << endl;

		cout << "=====================================" << endl << endl;
	
		// Check if C*C_inv gives identity
		I = Matrix_Product(C, C_inv);
		
		Print_Coeff(C, "C matrix");
		Print_Coeff(C_inv, "C_inv matrix");
		Print_Coeff(I, "I matrix");
		Print_Triangular(S, "S matrix");
		Print_Table(E, "Energies");
		
		I.clear();
	}
	
}

void Calculate_HKS () {

	vector< vector<double> > tmp;

	tmp = Matrix_Product(S,C);
	
	#pragma omp parallel for
	for (int i=0; i<n_mo; i++) {
		for (int j=0; j<n_mo; j++) {
			tmp[i][j]=tmp[i][j]*E[j];
		}
	}
	#pragma omp barrier
	
	H_KS = Matrix_Product(tmp,C_inv);

}

double Calculate_Transfer_Matrix (int mo_1, int mo_2) {

	double trans_mat[2][2], trans_mat_tmp[2][2], S12;

	//Sites Energies Calculations

	//E1
	trans_mat_tmp[0][0]= H_KS[mo_1][mo_1];

	//E2
	trans_mat_tmp[1][1]= H_KS[mo_2][mo_2];

	//Transfer Integrals

	//J12
	trans_mat_tmp[0][1]= H_KS[mo_1][mo_2];
	trans_mat_tmp[1][0]= H_KS[mo_2][mo_1];

	//S12
	S12= S[mo_1][mo_2];

	//printf("Transfer Matrix before othonormalization\n");
	//printf("%10.4lf   %10.4lf\n",trans_mat_tmp[0][0],trans_mat_tmp[0][1]);
	//printf("%10.4lf   %10.4lf\n",trans_mat_tmp[1][0],trans_mat_tmp[1][1]);
	//printf("S12=%lf\n",S12);
 
	//Effective parameters

	trans_mat[0][0]=(0.5 * (1.0/(1-S12*S12)) ) * (trans_mat_tmp[0][0] + trans_mat_tmp[1][1] - 2 * trans_mat_tmp[0][1] *S12 + (trans_mat_tmp[0][0] - trans_mat_tmp[1][1]) * sqrt (1 - S12*S12) );
	trans_mat[1][1]=(0.5 * (1.0/(1-S12*S12)) ) * (trans_mat_tmp[0][0] + trans_mat_tmp[1][1] - 2 * trans_mat_tmp[0][1] *S12 - (trans_mat_tmp[0][0] - trans_mat_tmp[1][1]) * sqrt (1 - S12*S12) );
	trans_mat[0][1]=(1.0/(1-S12*S12)) * (trans_mat_tmp[0][1] - 0.5 * (trans_mat_tmp[0][0] + trans_mat_tmp[1][1]) * S12);
	trans_mat[1][0]=(1.0/(1-S12*S12)) * (trans_mat_tmp[1][0] - 0.5 * (trans_mat_tmp[0][0] + trans_mat_tmp[1][1]) * S12);

	//printf("\nTransfer Matrix after othonormalization\n");
	//printf("%10.4lf   %10.4lf\n",trans_mat[0][0],trans_mat[0][1]);
	//printf("%10.4lf   %10.4lf\n",trans_mat[1][0],trans_mat[1][1]);
	
	// Return only J
	return(trans_mat[0][1]);

}

int main(int argc, char **argv) {

	time_t t_start, t_stop;
	string timeinfo_start, timeinfo_stop;
	
	// Get start time
	t_start = time(NULL);
	timeinfo_start = asctime(localtime(&t_start));
	
	//--------------------------------------------------
	
	string input_filename;
	int f, mol1, mol2;
	bool pytizi=false;

	int s;
	while ((s = getopt_long (argc, argv, "I:f:a:b:", NULL, NULL)) != -1){
		switch (s){
			case 'I':
				input_filename = optarg;
				break;
				
			case 'f':
				f = atoi(optarg);
				pytizi = true;
				break;
				
			case 'a':
				mol1 = atoi(optarg);
				break;
				
			case 'b':
				mol2 = atoi(optarg);
				break;
		}
	}

	Read_DIMER(input_filename, false);
	Calculate_HKS();
	
	double J_H, J_L;
	J_H = Calculate_Transfer_Matrix(homo_1, homo_2);
	J_L = Calculate_Transfer_Matrix(lumo_1, lumo_2);
	
	if (pytizi)
		cout << f << " " << mol1 << " " << mol2 << " " << scientific << setprecision(14) << J_H << " 1.0 1.0 " << J_L << " 1.0 1.0 " << endl;
		
	else
		cout << input_filename << " " << scientific << setprecision(14) << J_H << " " << J_L << endl;
	
	S.clear(); C.clear(); C_inv.clear(); I.clear(); E.clear();
	
	//--------------------------------------------------

	// Print final information
	t_stop = time(NULL);
	timeinfo_stop = asctime(localtime(&t_stop));
	
	//cout << "[INFO] Start time: " << timeinfo_start;	
	//cout << "[INFO] Stop time: " << timeinfo_stop;
	//cout << "[INFO] Calculation took " << t_stop-t_start << " seconds." << endl << endl;

	return(0);

}


