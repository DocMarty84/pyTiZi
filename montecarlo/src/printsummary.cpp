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

	OUT_TOT << output_folder << "/info_tot_" << charge.c_str() << "_" <<\
														F_dir.c_str() << ".dat";

	uF_x = F_x/F_norm;
	uF_y = F_y/F_norm;
	uF_z = F_z/F_norm;

	FILE * pFile;

	pFile = fopen(OUT_TOT.str().c_str(), "w");
	if (pFile==NULL) {
		int wait = 0; 
		while (wait<10 && pFile==NULL){
			cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) <<\
											" seconds to write a file" << endl;
			usleep(10*(wait+1)*(wait+1));
			pFile = fopen(OUT_TOT.str().c_str(), "w");
			wait++;
		}
		if (wait==10 && pFile==NULL){
			cerr << "[ERROR] Impossible to write " << OUT_TOT.str().c_str() <<\
														"! Exiting..." << endl;
			exit (1);
		}
	}
	fprintf(pFile,"Electric field unit vectors: (%f, %f, %f)\n\n",\
															uF_x, uF_y, uF_z);
	for (int i=0; i<n_frame; i++){
		fprintf(pFile,"Frame %d\n", i);
		for (int ii=0; ii<n_mol; ii++){
			fprintf(pFile,"Molecule %d | %d neighbors\n", mol_label[ii],\
												int(neigh_label[i][ii].size()));
			for (unsigned int jj=0; jj<neigh_label[i][ii].size(); jj++){
				fprintf(pFile,"%6d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %14.5"\
								"e\n", neigh_label[i][ii][jj],\
								d_x[i][ii][jj], d_y[i][ii][jj], d_z[i][ii][jj],\
								dE[i][ii][jj], J_H[i][ii][jj], J_L[i][ii][jj],\
								k[i][ii][jj]);
			}
		}
	}

	fclose(pFile);
}
