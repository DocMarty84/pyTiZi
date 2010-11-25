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

// Local libraries
#include "variables.h"

using namespace std;

// =============================================================================
// ------------------------ Coordinates transformations ------------------------
// =============================================================================

void Cartesian_To_Fractional(double* Dist_Cart, double* Dist_Frac, int i){
	double x = Dist_Cart[0];
	double y = Dist_Cart[1];
	double z = Dist_Cart[2];
	
	z = (z/temp_gamma_term[i]) / c[i];
	y = ((y-z*c[i]*temp_beta_term[i])/temp_gamma_sin[i]) / b[i];
	x = (x-y*b[i]*temp_gamma_cos[i]-z*c[i]*temp_beta_cos[i]) / a[i];
	Dist_Frac[0] = x; Dist_Frac[1] = y; Dist_Frac[2] = z;	
}

void Fractional_To_Cartesian(double* Dist_Frac, double* Dist_Cart, int i){
	double x = Dist_Frac[0];
	double y = Dist_Frac[1];
	double z = Dist_Frac[2];
		
	x = x*a[i] + y*b[i]*temp_gamma_cos[i] + z*c[i]*temp_beta_cos[i];
	y = y*b[i]*temp_gamma_sin[i] + z*c[i]*temp_beta_term[i];
	z = z*c[i]*temp_gamma_term[i];
	Dist_Cart[0] = x; Dist_Cart[1] = y; Dist_Cart[2] = z;
}
