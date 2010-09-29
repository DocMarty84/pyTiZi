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
vector< vector<double> > S, C, C_inv, I, E;

void Print_C(){

	int n_mo_rounded = (n_mo - (n_mo % 10))/10;
	
	for (int count_mo=0; count_mo<n_mo_rounded; count_mo++) {
		
		for(int i=0; i<n_mo; i++) {
			for(int j=0; j<11; j++) {

				if(j==1){
					printf("%6d%9.4lf", i+1, C[i][(j-1)+(count_mo*10)]);
				}
				
				else if(j>1){
					printf("%9.4lf", C[i][(j-1)+(count_mo*10)]);
				}
			}
			printf("\n");
		}
		printf("\n");
	}
	
	for(int i=0; i<n_mo; i++) {
		for(int j=0; j<(n_mo % 10)+1; j++) {
			
			if(j==1){
				printf("%6d%9.4lf", i+1, C[i][(j-1)+(n_mo_rounded*10)]);
			}
			
			else if(j>1){
				printf("%9.4f", C[i][(j-1)+(n_mo_rounded*10)]);
			}	
		}
		printf("\n");
	}
	printf("\n");
}

void Print_I(){

	int n_mo_rounded = (n_mo - (n_mo % 10))/10;
	
	for (int count_mo=0; count_mo<n_mo_rounded; count_mo++) {
		
		for(int i=0; i<n_mo; i++) {
			for(int j=0; j<11; j++) {

				if(j==1){
					printf("%6d%9.4lf", i+1, I[i][(j-1)+(count_mo*10)]);
				}
				
				else if(j>1){
					printf("%9.4lf", I[i][(j-1)+(count_mo*10)]);
				}
			}
			printf("\n");
		}
		printf("\n");
	}
	
	for(int i=0; i<n_mo; i++) {
		for(int j=0; j<(n_mo % 10)+1; j++) {
			
			if(j==1){
				printf("%6d%9.4lf", i+1, I[i][(j-1)+(n_mo_rounded*10)]);
			}
			
			else if(j>1){
				printf("%9.4f", I[i][(j-1)+(n_mo_rounded*10)]);
			}	
		}
		printf("\n");
	}
	printf("\n");
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

			if (str.find("Total nr. of (C)SFOs (summation over all irreps) :") != string::npos) {
				size_t length;
				char buffer[10];
				length = str.copy(buffer, 9, 51);
				buffer[length] = '\0';
				n_mo = atoi(buffer);
			}
			
			else if (str.find("MOs expanded in CFs+SFOs") != string::npos) {
				
				for(int i=0; i<n_mo; i++) {
					C.push_back( vector<double> ());
					C_inv.push_back( vector<double> ());
					I.push_back( vector<double> ());
					
					for(int j=0; j<n_mo; j++){
						C[i].push_back(0.0);
						C_inv[i].push_back(0.0);
						I[i].push_back(0.0);
					}
				}
				
				ap::real_2d_array a;
				a.setlength(n_mo, n_mo);
				
				int n_mo_rounded = (n_mo - (n_mo % 10))/10;
				
				for (int count_mo=0; count_mo<n_mo_rounded; count_mo++) {
					do {
						getline(input, str);
					} 
					while (str.find("CF+SFO") == string::npos);
					
					for(int i=0; i<n_mo; i++) {
						for(int j=0; j<11; j++) {
							input >> str;

							if(j>0){
								C[i][(j-1)+(count_mo*10)] = atof(str.c_str());
								a(i,(j-1)+(count_mo*10)) = atof(str.c_str());
							}		
						}
					}
				}
				
				do {
					getline(input, str);
				} 
				while (str.find("CF+SFO") == string::npos);
				
				for(int i=0; i<n_mo; i++) {
					for(int j=0; j<(n_mo % 10)+1; j++) {
						input >> str;

						if(j>0){
							C[i][(j-1)+(n_mo_rounded*10)] = atof(str.c_str());
							a(i,(j-1)+(n_mo_rounded*10)) = atof(str.c_str());
						}		
					}
				}
				
				if (print_result) {
					Print_C();
				}

				int info;
				matinvreport rep;
				rmatrixinverse(a, n_mo, info, rep);
				
				for(int i=0; i<n_mo; i++) {
					for(int j=0; j<n_mo; j++) {
						C_inv[i][j] = a(i,j);
					}
				}
				
				#pragma omp parallel for
				for (int i=0; i<n_mo; i++) {
					//cout << i << endl;
					double sum = 0.0;
					for (int j=0; j<n_mo; j++) {
						sum = 0.0;
						for (int k=0; k<n_mo; k++) {
							sum += C[i][k] * C_inv[k][j];
						}
						if (sum < 1e-12) {
							sum = 0.0;
						}
						I[i][j] = sum;
					}
				}
				#pragma omp barrier
				Print_I();
				
			}

		}
		input.close();
	}
	
	else{
		cerr << "Error opening " << input_file.c_str() << endl;
		exit(1);			
	}
}

/*
void calc_inv(bool print_results){

	ap::real_2d_array a;
	ap::real_2d_array b;
	ap::real_2d_array c;
	int info;
	matinvreport rep;

	double x;
	a.setlength(n, m);
	b.setlength(n, m);
	c.setlength(n, m);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
			x = rand();
			x = x/RAND_MAX;
			a(i,j) = int(x*100.0);
			b(i,j) = a(i,j);
			c(i,j) = 0.0;
		}
	}

	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		//cout << i << endl;
		double sum = 0.0;
		for (int j = 0; j < m; j++) {
			sum = 0.0;
			for (int k = 0; k < n; k++) {
				sum += a(i,k) * b(k,j);
			}
			if (sum < 1e-12) {
				sum = 0.0;
			}
			c(i,j) = sum;
		}
	}
	#pragma omp barrier

	rmatrixinverse(a, n, info, rep);

	if (print_results){
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < m; j++) {
				cout << b(i,j) << " ";
			}
			cout << endl;
		}
		cout << "----------------------" << endl;

		for(int i = 0; i < n; i++) {
		cout << "----------------------" << endl;
		
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < m; j++) {
				cout << c(i,j) << " ";
			}
			cout << endl;
		}
	}	for(int j = 0; j < m; j++) {
				cout << a(i,j) << " ";
			}
			cout << endl;
		}
		cout << "----------------------" << endl;
		
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < m; j++) {
				cout << c(i,j) << " ";
			}
			cout << endl;
		}
	}

}
*/


int main(int argc, char **argv) {

	string input_filename;

	int s;
	while ((s = getopt_long (argc, argv, "I:", NULL, NULL)) != -1){
		switch (s){
			case 'I':
				input_filename = optarg;
				break;
		}
	}

	Read_DIMER(input_filename, false);
	

}


