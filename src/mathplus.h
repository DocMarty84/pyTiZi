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

#ifndef _MATHPLUS_H
#define _MATHPLUS_H 1

// =============================================================================
// ----------------------------- Simple functions ------------------------------
// =============================================================================

// Return random number between 0 and 1
double Rand_0_1();

// Calculates n!
double Facto(int n);

// Product of two matrix
vector< vector<double> > Matrix_Product (vector< vector<double> > M1,\
                                         vector< vector<double> > M2,\
                                         unsigned int n_line_M1,\
                                         unsigned int n_line_M2);

// =============================================================================
// --------------------------------- Rotation ----------------------------------
// =============================================================================

// Calculates the rotation matrix related to any angle around any axis
// The rotation axis is defined by the vector rot_axis and the point origin
vector< vector<double> > Calcul_Rot_Matrix (double angle_deg,\
                                            vector<double> rot_axis,\
                                            vector<double> origin);

#endif // mathplus.h
