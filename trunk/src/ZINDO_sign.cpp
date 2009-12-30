#include <iostream> //Entrées-sorties standard
#include <fstream> //Entrées-sorties sur fichiers
#include <string> //Chaines de caracteres
//#include <limits> //Pour aller à la fin d'une ligne, par exemple
#include <iomanip> //Manipulation des sorties, par exemple format scientifique
#include <sstream> //Pour les conversion entre types de variables
#include <math.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void read_file(double *J_H, double *H_1, double *H_2, double *J_L, double *L_1, double *L_2, int n_frame, string input_file) {

	int tmp;
//	string input_file = "integrals.tmp"; 
	ifstream input(input_file.c_str(), ios::in);
	if(input)
	{
//		cout << "Transfer integrals file successfully opened! Now reading..." << endl;

		for(int i=0;i<n_frame;i++) input >> tmp >> J_H[i] >> H_1[i] >> H_2[i] >> J_L[i] >> L_1[i] >> L_2[i];
		input.close(); // Fermeture du fichier
//		cout << "Reading done!" << endl;
	}
	else
		cerr << "Error while opening trajectory file " << input_file.c_str() << endl;
}

void print_results(double *J_H, double *H_1, double *H_2, double *J_L, double *L_1, double *L_2, int n_frame) {
	for(int i=0;i<n_frame;i++){
	cout << i << " " << 1000*J_H[i]*(H_1[i]*H_2[i]/(fabs(H_1[i]*H_2[i]))) << " " << 1000*J_L[i]*(L_1[i]*L_2[i]/(fabs(L_1[i]*L_2[i]))) << endl;
	}
}

int main(int argc, char **argv) {

	int c;
	int n_frame = 0;
	string input_file;
	
  	while ((c = getopt_long (argc, argv, "f:n:", NULL, NULL)) != -1)
	{
      		switch (c)
		{
			case 'n':
				n_frame = atoi(optarg);
	  			break;
			case 'f':
				input_file = optarg;
	  			break;
		}
	}

//	cout << input_file << " " << n_frame << endl;
	if(n_frame == 0) {
		cout << "Frame number not specified: use -n option.\nAborting..." << endl;
		exit(1);
	}
	
	double *J_H, *H_1, *H_2, *J_L, *L_1, *L_2;
	
	J_H = new double[n_frame];
	H_1 = new double[n_frame];
	H_2 = new double[n_frame];
	J_L = new double[n_frame];
	L_1 = new double[n_frame];
	L_2 = new double[n_frame];
	
	read_file(J_H, H_1, H_2, J_L, L_1, L_2, n_frame, input_file);
	print_results(J_H, H_1, H_2, J_L, L_1, L_2, n_frame);
	
	delete[] J_H; delete[] H_1; delete[] H_2; 
	delete[] J_L; delete[] L_1; delete[] L_2;

	return 0;

}

