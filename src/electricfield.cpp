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
#include <vector>

// C libraries
#include <math.h>
#include <stdlib.h>

// Local libraries
#include "coordinates.h"
#include "variables.h"
#include "mathplus.h"

using namespace std;

// Calculates the electric field unit vector
void Calcul_F_vector(bool print_results) {
	
	if (F_dir.compare("a") == 0 || F_dir.compare("b") == 0 || F_dir.compare("c") == 0) {
		anisotropy = false;
		double *F_tmp_frac, *F_tmp_cart;
		F_tmp_frac = new double[3];
		F_tmp_cart = new double[3];

		if (F_dir.compare("a") == 0) {
			if (charge.compare("e") == 0) {
				F_tmp_frac[0] = 1.0;
			}
			else {
				F_tmp_frac[0] = -1.0;
			}
			F_tmp_frac[1] = 0.0;
			F_tmp_frac[2] = 0.0;
		}
			
		else if (F_dir.compare("b") == 0) {
			F_tmp_frac[0] = 0.0;
			if (charge.compare("e") == 0) {
				F_tmp_frac[1] = 1.0;
			}
			else {
				F_tmp_frac[1] = -1.0;
			}
			F_tmp_frac[2] = 0.0;
		}
			
		else if (F_dir.compare("c") == 0) {
			F_tmp_frac[0] = 0.0;
			F_tmp_frac[1] = 0.0;
			if (charge.compare("e") == 0) {
				F_tmp_frac[2] = 1.0;
			}
			else {
				F_tmp_frac[2] = -1.0;
			}
		}
					
		Fractional_To_Cartesian(F_tmp_frac, F_tmp_cart, 0);
		
		double F_norm_tmp = sqrt(pow(F_tmp_cart[0],2) + pow(F_tmp_cart[1],2) + pow(F_tmp_cart[2],2));
		F_x_list.push_back(F_tmp_cart[0]/F_norm_tmp);
		F_y_list.push_back(F_tmp_cart[1]/F_norm_tmp);
		F_z_list.push_back(F_tmp_cart[2]/F_norm_tmp);
		F_angle_list.push_back(0.0);
		
		delete [] F_tmp_frac; delete [] F_tmp_cart;
		
	}
	
	else {
		anisotropy = true;
		double *v1, *v2, *v3, *v_tmp_frac, *v_tmp_cart;
		v1 = new double[3];
		v2 = new double[3];
		v3 = new double[3];
		v_tmp_frac = new double[3];
		v_tmp_cart = new double[3];
	
		vector< vector<double> > Rot_Matrix (4, vector<double> (4, 0.0));
		vector< vector<double> > F_0_tmp_cart (4, vector<double> (1, 0.0));
		vector< vector<double> > F_tmp_cart (4, vector<double> (1, 0.0));
		
		if (F_dir.compare("ab") == 0) {
			
			// Get two axis parallel to a and b
			v_tmp_frac[0] = 1.0;
			v_tmp_frac[1] = 0.0;
			v_tmp_frac[2] = 0.0;
			Fractional_To_Cartesian(v_tmp_frac, v_tmp_cart, 0);
			v1[0] = v_tmp_cart[0];
			v1[1] = v_tmp_cart[1];
			v1[2] = v_tmp_cart[2];
			
			v_tmp_frac[0] = 0.0;
			v_tmp_frac[1] = 1.0;
			v_tmp_frac[2] = 0.0;
			Fractional_To_Cartesian(v_tmp_frac, v_tmp_cart, 0);
			v2[0] = v_tmp_cart[0];
			v2[1] = v_tmp_cart[1];
			v2[2] = v_tmp_cart[2];
		}
		
		else if (F_dir.compare("ac") == 0) {
			
			// Get two axis parallel to a and c
			v_tmp_frac[0] = 1.0;
			v_tmp_frac[1] = 0.0;
			v_tmp_frac[2] = 0.0;
			Fractional_To_Cartesian(v_tmp_frac, v_tmp_cart, 0);
			v1[0] = v_tmp_cart[0];
			v1[1] = v_tmp_cart[1];
			v1[2] = v_tmp_cart[2];
			
			v_tmp_frac[0] = 0.0;
			v_tmp_frac[1] = 0.0;
			v_tmp_frac[2] = 1.0;
			Fractional_To_Cartesian(v_tmp_frac, v_tmp_cart, 0);
			v2[0] = v_tmp_cart[0];
			v2[1] = v_tmp_cart[1];
			v2[2] = v_tmp_cart[2];
		}

		else if (F_dir.compare("bc") == 0) {
			
			// Get two axis parallel to b and c
			v_tmp_frac[0] = 0.0;
			v_tmp_frac[1] = 1.0;
			v_tmp_frac[2] = 0.0;
			Fractional_To_Cartesian(v_tmp_frac, v_tmp_cart, 0);
			v1[0] = v_tmp_cart[0];
			v1[1] = v_tmp_cart[1];
			v1[2] = v_tmp_cart[2];
			
			v_tmp_frac[0] = 0.0;
			v_tmp_frac[1] = 0.0;
			v_tmp_frac[2] = 1.0;
			Fractional_To_Cartesian(v_tmp_frac, v_tmp_cart, 0);
			v2[0] = v_tmp_cart[0];
			v2[1] = v_tmp_cart[1];
			v2[2] = v_tmp_cart[2];
		}
		
		// Get the third axis perpendicular to the plane
		v3[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
		v3[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
		v3[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);

		for (int angle=0; angle<360; angle = angle + 15) {
			// Get the rotation matrix
			double theta_deg = float(angle);
			vector<double> origin (3, 0.0);
			vector<double> rot_axis (3, 0.0);
			rot_axis[0] = v3[0];
			rot_axis[1] = v3[1];
			rot_axis[2] = v3[2];
			Rot_Matrix = Calcul_Rot_Matrix(theta_deg, rot_axis, origin);
			
			// Get the rotated electric field axis
			F_0_tmp_cart[0][0] = v1[0];
			F_0_tmp_cart[1][0] = v1[1];
			F_0_tmp_cart[2][0] = v1[2];
			F_0_tmp_cart[3][0] = 1.0;
			F_tmp_cart = Matrix_Product(Rot_Matrix, F_0_tmp_cart, 4, 4);
			double F_norm_tmp = sqrt(pow(F_0_tmp_cart[0][0],2) + pow(F_0_tmp_cart[1][0],2) + pow(F_0_tmp_cart[2][0],2));
			
			F_angle_list.push_back(theta_deg);
			F_x_list.push_back(F_tmp_cart[0][0]/F_norm_tmp);
			F_y_list.push_back(F_tmp_cart[1][0]/F_norm_tmp);
			F_z_list.push_back(F_tmp_cart[2][0]/F_norm_tmp);
		}	
		
		delete [] v1; delete [] v2; delete [] v3;
		delete [] v_tmp_frac; delete [] v_tmp_cart;
	}
	
	for (unsigned int m=0; m<F_angle_list.size(); m++){
		if (fabs(F_x_list[m]) < 1E-10) F_x_list[m] = 0.0;
		if (fabs(F_y_list[m]) < 1E-10) F_y_list[m] = 0.0;
		if (fabs(F_z_list[m]) < 1E-10) F_z_list[m] = 0.0;
		
		F_x_list[m] = F_x_list[m] * F_norm; 
		F_y_list[m] = F_y_list[m] * F_norm; 
		F_z_list[m] = F_z_list[m] * F_norm;
	}
	
	if (print_results) {
		for (unsigned int m=0; m<F_angle_list.size(); m++){
			cout << F_dir << " plane: angle = " << F_angle_list[m] << " | F_x = " << F_x_list[m] << " | F_y = " << F_y_list[m] << " | F_z = " << F_z_list[m] << endl;
		}
	}
}
