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

#ifndef _DISPATCH_H
#define _DISPATCH_H 1

// Choose a starting molecule randomly
int Choose_Mol_RND(int frame);

// ===========================================================================================================
// ------------------------------------------- Dispatch Charges ----------------------------------------------
// ===========================================================================================================

// Dispatch the charges randomly in the grid
void Dispatch_Mol_RND(int frame, vector< vector<bool> > grid_occ, int *pos);

// Dispatch the charges randomly in the grid, but in a specific layer
void Dispatch_Mol_RND_layer(int frame, vector< vector<bool> > grid_occ, int *pos);

// Dispatch the charges randomly in the fist "column" of the grid
void Dispatch_Mol_begin(int frame, vector< vector<bool> > grid_occ, int *pos);

// Dispatch the charges randomly in the fist "column" of the grid, 
// but in a specific layer
void Dispatch_Mol_begin_layer(int frame, vector< vector<bool> > grid_occ, int *pos);

#endif // dispatch.h
