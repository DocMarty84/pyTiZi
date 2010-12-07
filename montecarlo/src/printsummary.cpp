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
#include <iostream> 
#include <fstream> 
#include <string> 
#include <sstream> 
#include <vector>

// C libraries
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Local libraries
#include "constants.h"
#include "variables.h"

using namespace std;


// ===========================================================================================================
// --------------------------------------------- Print part --------------------------------------------------
// ===========================================================================================================

// Print all the informations used to calculate k (before running any simulation)
void Print_Summary_Beginning(string output_folder) {
	stringstream OUT_TOT;
	double uF_x, uF_y, uF_z;

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
				fprintf(pFile,"%6d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %14.5e\n", neigh_label[i][ii][jj],\
												d_x[i][ii][jj], d_y[i][ii][jj], d_z[i][ii][jj],\
												dE[i][ii][jj], J_H[i][ii][jj], J_L[i][ii][jj], k[i][ii][jj]);
			}
		}
	}

	fclose(pFile);
}

// Print some information for each try of the MC simulation
void Print_Summary_Try(string output_folder, int i, int charge_try, double total_dist_try,\
																					double total_time_try) {
	stringstream OUT_SIMU_FRAME;
	double uF_x, uF_y, uF_z;

	OUT_SIMU_FRAME << output_folder << "/simu_" << charge.c_str() << "_" << F_dir.c_str() << ".out";

	uF_x = F_x/F_norm;
	uF_y = F_y/F_norm;
	uF_z = F_z/F_norm;

	FILE * pFile;
	
	pFile=fopen(OUT_SIMU_FRAME.str().c_str(), "a");
	if (pFile==NULL) {
		int wait = 0; 
		while (wait<10 && pFile==NULL){
			cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) << " seconds to write a file" << endl;
			usleep(10*(wait+1)*(wait+1));
			pFile=fopen(OUT_SIMU_FRAME.str().c_str(), "a");
			wait++;
		}
		if (wait==10 && pFile==NULL){
			cerr << "[ERROR] Impossible to write " << OUT_SIMU_FRAME.str().c_str() << "! Exiting..." << endl;
			exit (1);
		}
	}
	fprintf(pFile,"-------------------------------------------------------------------------------\n");
	fprintf(pFile,"Frame = %d\n", i);
	fprintf(pFile,"Electric_Field_Angle = %d\n", int(F_angle));
	fprintf(pFile,"Electric_Field_Unit_Vectors = (%f, %f, %f)\n", uF_x, uF_y, uF_z);
	fprintf(pFile,"Number_of_Charges = %d\n", n_charges);
	fprintf(pFile,"Density_of_Charges = %.5e charges/cm3\n", double(n_charges)/(vol_box[i]*n_mini_grid_a*n_mini_grid_b*n_mini_grid_c*1e-24));
	fprintf(pFile,"Time_try_%d = %e\n", (charge_try/n_charges), total_time_try/(charge_try/n_charges));
	fprintf(pFile,"Distance_try_%d = %e\n", (charge_try/n_charges), total_dist_try/(charge_try/n_charges));
	fprintf(pFile,"Mu_try_%d = %lf\n", (charge_try/n_charges), total_dist_try/(total_time_try*F_norm));
	fclose(pFile);
}
