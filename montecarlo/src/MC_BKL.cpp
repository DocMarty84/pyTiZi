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

#define _INSIDE_MC_BKL

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
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <omp.h>

// Local libraries
#include "constants.h"
#include "variables.h"
#include "clear.h"
#include "coordinates.h"
#include "mathplus.h"
#include "read.h"
#include "buildgrid.h"
#include "electricfield.h"
#include "distance.h"
#include "deltae.h"
#include "deltav.h"
#include "transferrate.h"
#include "dispatch.h"
#include "printsummary.h"
#include "mcbkl.h"
#include "mcfrm.h"

using namespace std;

// =============================================================================
// ----------------------- Monte-Carlo related functions -----------------------
// =============================================================================

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
					Clear_All();
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
		Clear_All();
		exit(1);
	}
	
	if (F_dir.compare("a") == 0 || F_dir.compare("b") == 0 || F_dir.compare("c") == 0)
		cout << "[INFO] Electric field is along the '" << F_dir << "' direction." << endl;
		
	else if (F_dir.compare("ab") == 0 || F_dir.compare("ac") == 0 || F_dir.compare("bc") == 0)
		cout << "[INFO] Electric field will probe the anisotropy in the '" << F_dir << "'plane." << endl;
		
	else {
		cerr << "[ERROR] Direction of the electric field not specified! Please use -d {a,b,c,ab,ac,bc}. Exiting..." << endl;
		Clear_All();
		exit(1);
	}

	// Read the required files
	Read_MC(input_file, input_folder, false);
	Read_CELL(input_file, input_folder, false);
	Read_CM(input_file, input_folder, false);
	Read_E_av(input_file, input_folder, false);

	// Set up lambda_i
	if (charge.compare("e") == 0) {
		LAMBDA_I = LAMBDA_I_E;
	}
	else {
		LAMBDA_I = LAMBDA_I_H;
	}

	// Calculates Electric field vectors
	Calcul_F_vector(false);
	
	// Calculates distances and DeltaE
	Calcul_Dist(false);
	Calcul_DeltaE(false);

	// Build the grid
	Build_Grid(true);
	
	// Save the triangular matrix
	vector< vector< vector<int> > > neigh_label_ref = neigh_label;
	vector< vector< vector<int> > > neigh_jump_vec_a_ref = neigh_jump_vec_a, neigh_jump_vec_b_ref = neigh_jump_vec_b, neigh_jump_vec_c_ref = neigh_jump_vec_c;
	vector< vector< vector<double> > > J_H_ref = J_H, J_L_ref = J_L;
	vector< vector< vector<double> > > d_x_ref = d_x, d_y_ref = d_y, d_z_ref = d_z;
	vector< vector< vector<double> > > dE_ref = dE;				
	
	//#pragma omp parallel for private(F_x, F_y, F_z, F_angle, F_x_list, F_y_list, F_z_list, F_angle_list, k, k_inv, grid_occ) //if (n_frame < int(omp_get_max_threads()) && anisotropy == true)
	for (unsigned int m=0; m<F_angle_list.size(); m++) {
		
		F_x = F_x_list[m];
		F_y = F_y_list[m];
		F_z = F_z_list[m];
		F_angle = F_angle_list[m];

		// Calculate transfer rates and the full matrix, mostly for information
		Marcus_Levich_Jortner_CST(); // Calculates constants
		Calcul_k(false);
		Full_Matrix();

		// Print a summary
		Print_Summary(output_folder);
		
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
	}
	//#pragma omp barrier

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
