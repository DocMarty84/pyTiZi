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
#include <numeric> 
#include <string> 
#include <sstream> 
#include <vector>

// C libraries
#include <math.h>
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
											dE_box[i][ii][jj], J_H[i][ii][jj], J_L[i][ii][jj], k[i][ii][jj]);
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

	OUT_SIMU_FRAME << output_folder << "/simu_" << charge.c_str() << "_" << F_dir.c_str()  << "_f_" << i <<\
																									".out";

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
	fprintf(pFile,"Frame           = %d\n", i);
	fprintf(pFile,"Time_try_%d     = %e\n", (charge_try/n_charges), total_time_try/(charge_try/n_charges));
	fprintf(pFile,"Distance_try_%d = %e\n", (charge_try/n_charges), total_dist_try/(charge_try/n_charges));
	fprintf(pFile,"Mu_try_%d       = %lf\n", (charge_try/n_charges), total_dist_try/(total_time_try*F_norm));
	fclose(pFile);
}

// Print some information at the end of each frame calculation
void Print_Summary_Frame(string output_folder, int i, double total_dist_try, double total_time_try,\
																				vector <double> mu_frame) {
	stringstream OUT_SIMU_FRAME;
	double uF_x, uF_y, uF_z;

	OUT_SIMU_FRAME << output_folder << "/simu_" << charge.c_str() << "_" << F_dir.c_str()  << "_f_" << i <<\
																									".out";

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
	fprintf(pFile,"\n\n###############################################################################\n");
	fprintf(pFile,"##                                                                           ##\n");
	fprintf(pFile,"##                          Summary for the Frame%6d                      ##\n", i);
	fprintf(pFile,"##                                                                           ##\n");
	fprintf(pFile,"###############################################################################\n\n");
	fprintf(pFile,"Global Parameters\n");
	fprintf(pFile,"-----------------\n");
	fprintf(pFile,"Frame = %d\n", i);
	fprintf(pFile,"a     = %lf\n", a[i]);
	fprintf(pFile,"b     = %lf\n", b[i]);
	fprintf(pFile,"c     = %lf\n", c[i]);
	fprintf(pFile,"alpha = %lf\n", alpha_deg[i]);
	fprintf(pFile,"beta  = %lf\n", beta_deg[i]);
	fprintf(pFile,"gamma = %lf\n", gamma_deg[i]);
	fprintf(pFile,"volume unit cell = %e Angstrom3\n\n", vol_box[i]);
	
	fprintf(pFile,"Monte-Carlo Parameters\n");
	fprintf(pFile,"----------------------\n");	
	fprintf(pFile,"Number of tries             = %d\n", n_try);
	fprintf(pFile,"Number of boxes along a     = %d\n", n_box_a);
	fprintf(pFile,"Number of boxes along b     = %d\n", n_box_b);
	fprintf(pFile,"Number of boxes along c     = %d\n", n_box_c);
	fprintf(pFile,"Volume of grid              = %e Angstrom3\n", vol_box[i]*n_box);
	fprintf(pFile,"Electric Field Norm         = %e V/cm\n", F_norm);
	fprintf(pFile,"Electric Field Direction    = %s \n", F_dir.c_str());
	fprintf(pFile,"Electric Field Angle        = %d degrees\n", int(F_angle));
	fprintf(pFile,"Electric Field Unit Vectors = (%f, %f, %f)\n", uF_x, uF_y, uF_z);
	fprintf(pFile,"Type of Charges             = %s\n", charge.c_str());
	fprintf(pFile,"Number of Charges           = %d\n", n_charges);
	fprintf(pFile,"Number of Charges per Site  = %e charges/site\n", double(n_charges)/(n_mol*n_box));
	fprintf(pFile,"Density of Charges          = %e charges/cm3\n\n", double(n_charges)/(vol_box[i]*n_box*\
																									1e-24));
	fprintf(pFile,"Monte-Carlo Results\n");
	fprintf(pFile,"-------------------\n");																						
	fprintf(pFile,"Average Time     = %e seconds\n", total_time_try/double(n_try*n_charges));
	fprintf(pFile,"Average Distance = %e centimeters\n", total_dist_try/double(n_try*n_charges));
	fprintf(pFile,"Mobility         = %e cm2/Vs\n\n", mu_frame[i]);

	// Print the properties of each charge
	if (VERB_RECORD) {
		vector< vector< double> > chrg_E_electrostatic_sq(chrg_E_electrostatic[i]);
		vector< vector< double> > chrg_E_0_sq(chrg_E_0[i]); 
		vector< vector< double> > chrg_E_1_sq(chrg_E_1[i]);
		vector<double>::iterator chrg_E_electrostatic_it, chrg_E_0_it, chrg_E_1_it;
		
		vector <double> chrg_E_electrostatic_av, chrg_E_0_av, chrg_E_1_av;
		vector <double> chrg_E_electrostatic_sq_av, chrg_E_0_sq_av, chrg_E_1_sq_av;
		vector <double> chrg_E_electrostatic_var, chrg_E_0_var, chrg_E_1_var;
		
		for (unsigned int charge_i=0; charge_i<chrg_E_electrostatic.size(); charge_i++){
			// Calculate square of each value (should work for huge tables, I hope)
			for (chrg_E_electrostatic_it=chrg_E_electrostatic_sq[charge_i].begin(); chrg_E_electrostatic_it <\
											chrg_E_electrostatic_sq[charge_i].end(); chrg_E_electrostatic_it++) {
				*chrg_E_electrostatic_it = pow(*chrg_E_electrostatic_it, 2);
			}
			for (chrg_E_0_it=chrg_E_0_sq[charge_i].begin(); chrg_E_0_it < chrg_E_0_sq[charge_i].end();\
																								chrg_E_0_it++) {
				*chrg_E_0_it = pow(*chrg_E_0_it, 2);
			}
			for (chrg_E_1_it=chrg_E_1_sq[charge_i].begin(); chrg_E_1_it < chrg_E_1_sq[charge_i].end();\
																								chrg_E_1_it++) {
				*chrg_E_1_it = pow(*chrg_E_1_it, 2);
			}
			
			// Calculate average
			chrg_E_electrostatic_av.push_back(accumulate(chrg_E_electrostatic[i][charge_i].begin(),\
						chrg_E_electrostatic[i][charge_i].end(), 0.0)/chrg_E_electrostatic[i][charge_i].size());
			chrg_E_0_av.push_back(accumulate(chrg_E_0[i][charge_i].begin(),\
												chrg_E_0[i][charge_i].end(), 0.0)/chrg_E_0[i][charge_i].size());
			chrg_E_1_av.push_back(accumulate(chrg_E_1[i][charge_i].begin(),\
												chrg_E_1[i][charge_i].end(), 0.0)/chrg_E_1[i][charge_i].size());
														
			// Calculate average of the squared values
			chrg_E_electrostatic_sq_av.push_back(accumulate(chrg_E_electrostatic_sq[charge_i].begin(),\
						chrg_E_electrostatic_sq[charge_i].end(), 0.0)/chrg_E_electrostatic_sq[charge_i].size());
			chrg_E_0_sq_av.push_back(accumulate(chrg_E_0_sq[charge_i].begin(),\
												chrg_E_0_sq[charge_i].end(), 0.0)/chrg_E_0_sq[charge_i].size());
			chrg_E_1_sq_av.push_back(accumulate(chrg_E_1_sq[charge_i].begin(),\
												chrg_E_1_sq[charge_i].end(), 0.0)/chrg_E_1_sq[charge_i].size());

			// Calculate the variance									
			chrg_E_electrostatic_var.push_back(chrg_E_electrostatic_sq_av.back() -\
																		pow(chrg_E_electrostatic_av.back(), 2));
			chrg_E_0_var.push_back(chrg_E_0_sq_av.back() - pow(chrg_E_0_av.back(), 2));
			chrg_E_1_var.push_back(chrg_E_1_sq_av.back() - pow(chrg_E_1_av.back(), 2));
		}

		fprintf(pFile,"Charge properties\n");
		fprintf(pFile,"-----------------\n");
		fprintf(pFile,"Charge_number V_average V_variance E_0_average E_0_variance E_1_average E_1_variance "\
																		"Av_Travel_time Av_Distance Mobility\n");	
		for (unsigned int charge_i=0; charge_i<chrg_E_electrostatic[i].size(); charge_i++){
			fprintf(pFile,"%d %e %e %e %e %e %e %e %e %e\n", charge_i, chrg_E_electrostatic_av[charge_i],\
						chrg_E_electrostatic_var[charge_i], chrg_E_0_av[charge_i], chrg_E_0_var[charge_i],\
						chrg_E_1_av[charge_i], chrg_E_1_var[charge_i],\
						chrg_total_time[i][charge_i]/double(n_try), chrg_total_dist[i][charge_i]/double(n_try),\
						chrg_total_dist[i][charge_i]/(chrg_total_time[i][charge_i]*F_norm));
		}
	}
	
	else {
		fprintf(pFile,"Charge properties\n");
		fprintf(pFile,"-----------------\n");
		fprintf(pFile,"Charge_number Av_Travel_time Av_Distance Mobility\n");	
		for (unsigned int charge_i=0; charge_i<chrg_E_electrostatic[i].size(); charge_i++){
			fprintf(pFile,"%d %e %e %e\n", charge_i, chrg_total_time[i][charge_i]/double(n_try),\
										chrg_total_dist[i][charge_i]/double(n_try),\
										chrg_total_dist[i][charge_i]/(chrg_total_time[i][charge_i]*F_norm));
		}	
	}

	// Print the properties of each site of the grid
	fprintf(pFile,"\nGrid properties\n");
	fprintf(pFile,"---------------\n");	

	fprintf(pFile,"Box Molecule   X Y Z   E_0 E_1 Probability\n");
	fprintf(pFile," x      x    Angstrom     eV         x \n");

	double grid_probability_norm = 0.0;
	for (int x=0; x<n_box; x++){
		for (int ii=0; ii<n_mol; ii++){
			for (unsigned int charge_i=0; charge_i < grid_probability[i][x][ii].size(); charge_i++) {
				grid_probability_norm += grid_probability[i][x][ii][charge_i];
			}
		}
	}
	
	for (int x=0; x<n_box; x++){
		for (int ii=0; ii<n_mol; ii++){
			
			double grid_probability_av = 0.0;
			for (unsigned int charge_i=0; charge_i < grid_probability[i][x][ii].size(); charge_i++) {
				grid_probability_av += grid_probability[i][x][ii][charge_i];
			}
			grid_probability_av /= (double(grid_probability[i][x][ii].size())*grid_probability_norm);
			
			fprintf(pFile, "%d %d %lf %lf %lf %e %e %e\n", x, ii,\
												grid_x[i][x][ii], grid_y[i][x][ii], grid_z[i][x][ii],\
												grid_E_0[i][x][ii], grid_E_1[i][x][ii], grid_probability_av);
		}
	}
	fclose(pFile);
}

// Print some information at the end of the simulation
void Print_Summary_Final(string output_folder, double mu_moy) {
	stringstream OUT_SIMU_FINAL;
	double uF_x, uF_y, uF_z;

	OUT_SIMU_FINAL << output_folder << "/simu_" << charge.c_str() << "_" << F_dir.c_str() << ".out";

	uF_x = F_x/F_norm;
	uF_y = F_y/F_norm;
	uF_z = F_z/F_norm;

	FILE * pFile;
		
	pFile=fopen(OUT_SIMU_FINAL.str().c_str(), "a");
	if (pFile==NULL) {
		int wait = 0; 
		while (wait<10 && pFile==NULL){
			cerr << "[ERROR] Waiting " << 10*(wait+1)*(wait+1) << " seconds to write a file" << endl;
			usleep(10*(wait+1)*(wait+1));
			pFile=fopen(OUT_SIMU_FINAL.str().c_str(), "a");
			wait++;
		}
		if (wait==10 && pFile==NULL){
			cerr << "[ERROR] Impossible to write " << OUT_SIMU_FINAL.str().c_str() << "! Exiting..." << endl;
			exit (1);
		}
	}
	fprintf(pFile,"\n###############################################################################\n");
	fprintf(pFile,"##                                                                           ##\n");
	fprintf(pFile,"##                                 Final Summary                             ##\n");
	fprintf(pFile,"##                                                                           ##\n");
	fprintf(pFile,"###############################################################################\n\n");
	fprintf(pFile,"Global Parameters\n");
	fprintf(pFile,"-----------------\n");
	fprintf(pFile,"Number of Frames            = %d\n", n_frame);
	fprintf(pFile,"Number of Tries             = %d\n", n_try);
	fprintf(pFile,"Electric Field Angle        = %d\n", int(F_angle));
	fprintf(pFile,"Electric Field Unit Vectors = (%f, %f, %f)\n", uF_x, uF_y, uF_z);
	fprintf(pFile,"Number of Charges           = %d\n", n_charges);
	fprintf(pFile,"Number of Charges per Site  = %e charges/site\n", double(n_charges)/(n_mol*n_box));
	fprintf(pFile,"Density of Charges          = %e charges/cm3\n", double(n_charges)/(vol_box[0]*\
																								n_box*1e-24));
	fprintf(pFile,"Average Mobility            = %e cm2/Vs\n", mu_moy);
	fclose(pFile);
	

}
