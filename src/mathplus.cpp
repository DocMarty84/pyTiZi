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
#include "constants.h"
#include "clear.h"

using namespace std;

// =============================================================================
// ----------------------------- Simple functions ------------------------------
// =============================================================================

// Return random number between 0 and 1
double Rand_0_1(){
	double x;
	do {
		x = rand();
		x = x/RAND_MAX;
	}
	while(x==0);
	
	return x;
}

// Calculates n!
double Facto(int n) {
	if (n < 2) 
		return 1.0;
  
  	else {
		double x = 1.0;
		for (int i = 2; i <= n; i++) {
			x *= i;
		}
		return x;
	}
}

// Product of two matrix
vector< vector<double> > Matrix_Product (vector< vector<double> > M1,\
                                         vector< vector<double> > M2,\
                                         unsigned int n_line_M1,\
                                         unsigned int n_line_M2) {
	if (M1[0].size() != n_line_M2) {
		cerr << "[ERROR] Could not multiply matrix!\nExiting..." << endl;
		Clear_All();
		exit(1);
	}

	vector< vector<double> > M_res;
	
	for(unsigned int i=0; i<n_line_M1; i++) {
		M_res.push_back( vector<double> ());
		
		for(unsigned int j=0; j<M2[0].size(); j++){
			M_res[i].push_back(0.0);
		}
	}

	for (unsigned int i=0; i<n_line_M1; i++) {
		double sum = 0.0;
		for (unsigned int j=0; j<M2[0].size(); j++) {
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
	
	return(M_res);

}

// =============================================================================
// --------------------------------- Rotation ----------------------------------
// =============================================================================

// Calculates the rotation matrix related to any angle around any axis
// The rotation axis is defined by the vector rot_axis and the point origin
vector< vector<double> > Calcul_Rot_Matrix (double angle_deg,\
                                            vector<double> rot_axis,\
                                            vector<double> origin) {
	
	vector< vector<double> > Rot_Matrix (4, vector<double> (4, 0.0));
	
	double u, v, w, u2, v2, w2, L, L2;
	double i,j,k;
	double angle_rad;
	
	u = rot_axis[0];
	v = rot_axis[1];
	w = rot_axis[2];
	u2 = u*u;
	v2 = v*v;
	w2 = w*w;
	L2 = u2 + v2 + w2;
	L = sqrt(L2);
	
	i = origin[0];
	j = origin[1];
	k = origin[2];
	
	angle_rad = (PI/180.0) * angle_deg;
	
	Rot_Matrix[0][0] = (u2 + (v2 + w2) * cos(angle_rad))/L2;
	Rot_Matrix[0][1] = (u * v * (1.0 - cos(angle_rad)) -\
	                    w * L * sin(angle_rad))/L2;
	Rot_Matrix[0][2] = (u * w * (1.0 - cos(angle_rad)) +\
	                    v * L * sin(angle_rad))/L2;
	Rot_Matrix[0][3] = (i*(v2 + w2) - u*(j*v + k*w) + (u*(j*v + k*w) -\
	         i*(v2+w2)) * cos(angle_rad) + (j*w - k*v) * L * sin(angle_rad))/L2;
	Rot_Matrix[1][0] = (u * v * (1.0 - cos(angle_rad)) +\
	                    w * L * sin(angle_rad))/L2;
	Rot_Matrix[1][1] = (v2 + (u2 + w2) * cos(angle_rad))/L2;
	Rot_Matrix[1][2] = (v * w * (1.0 - cos(angle_rad)) -\
	                    u * L * sin(angle_rad))/L2;
	Rot_Matrix[1][3] = (j*(u2 + w2) - v*(i*u + k*w) + (v*(i*u + k*w) -\
	         j*(u2+w2)) * cos(angle_rad) + (k*u - i*w) * L * sin(angle_rad))/L2;
	Rot_Matrix[2][0] = (u * w * (1.0 - cos(angle_rad)) -\
	                    v * L * sin(angle_rad))/L2;
	Rot_Matrix[2][1] = (v * w * (1.0 - cos(angle_rad)) +\
	                    u * L * sin(angle_rad))/L2;
	Rot_Matrix[2][2] = (w2 + (u2 + v2) * cos(angle_rad))/L2;
	Rot_Matrix[2][3] = (k*(u2 + v2) - w*(i*u + j*v) + (w*(i*u + j*v) -\
	         k*(u2+v2)) * cos(angle_rad) + (i*v - j*u) * L * sin(angle_rad))/L2;
	Rot_Matrix[3][0] = 0.0;
	Rot_Matrix[3][1] = 0.0;
	Rot_Matrix[3][2] = 0.0;
	Rot_Matrix[3][3] = 1.0;
	
	return Rot_Matrix ;
	
}
