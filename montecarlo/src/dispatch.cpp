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
#include <vector>

// C libraries
#include <stdlib.h>

// Local libraries
#include "variables.h"

using namespace std;

// Choose a starting molecule randomly
int Choose_Mol_RND(int frame){ 
	int mol; double k_sum;
	
	do {
		mol = rand()%n_mol;
		k_sum = 0.0;
		
		for(unsigned int jj=0; jj<neigh_label[frame][mol].size(); jj++) {
			k_sum = k_sum + k_inv[frame][mol][jj];
		}
	}
	while(neigh_label[frame][mol].size()==0 || neigh_label[frame][mol].size() == 0);

	return mol;
}

// ===========================================================================================================
// ------------------------------------------- Dispatch Charges ----------------------------------------------
// ===========================================================================================================

void Dispatch_Mol_RND(int frame, vector< vector<bool> > grid_occ, int *pos){ 
	
	int pos_a=0, pos_b=0, pos_c=0, mol=0, box=0; 
	
	//double k_sum;
	
	do {
		pos_a = rand()%n_box_a;
		pos_b = rand()%n_box_b;
		pos_c = rand()%n_box_c;

		mol = rand()%n_mol;
		
		//k_sum = 0.0;
		//
		//for(unsigned int jj=0; jj<neigh_label[frame][mol].size(); jj++) {
		//	k_sum = k_sum + k_inv[frame][mol][jj];
		//}
		
		// Find the number of the corresponding mini-grid
		for (unsigned int x=0; x<box_neigh_label.size(); x++){
			if (pos_a == box_a[x] && pos_b == box_b[x] && pos_c == box_c[x]){
				box = x;
				break;
			}
		}
	}
	while(neigh_label[frame][mol].size() == 0 || grid_occ[box][mol] == true);
	
	pos[0] = mol;
	pos[1] = box;
}

void Dispatch_Mol_RND_layer(int frame, vector< vector<bool> > grid_occ,	int *pos){ 
	
	int pos_a=0, pos_b=0, pos_c=0, mol=0, box=0; 
	
	//double k_sum;
	
	do {
		pos_a = rand()%n_box_a;
		pos_b = rand()%n_box_b;
		pos_c = rand()%n_box_c;

		//mol = rand()%n_mol;
		mol = mol_layer[frame][layer][rand()%mol_layer[frame][layer].size()];
	
		//k_sum = 0.0;
		//
		//for(unsigned int jj=0; jj<neigh_label[frame][mol].size(); jj++) {
		//	k_sum = k_sum + k_inv[frame][mol][jj];
		//}
		
		// Find the number of the corresponding mini-grid
		for (unsigned int x=0; x<box_neigh_label.size(); x++){
			if (pos_a == box_a[x] && pos_b == box_b[x] && pos_c == box_c[x]){
				box = x;
				break;
			}
		}
	}
	while(neigh_label[frame][mol].size() == 0 || grid_occ[box][mol] == true);
	
	pos[0] = mol;
	pos[1] = box;
}

void Dispatch_Mol_begin(int frame, vector< vector<bool> > grid_occ, int *pos){ 
	
	int pos_a=0, pos_b=0, pos_c=0, mol=0, box=0; 
	
	//double k_sum;
	
	do {
		if (F_dir.compare("a") == 0) {
			pos_a = 0;
			pos_b = rand()%n_box_b;
			pos_c = rand()%n_box_c;
		}
		
		else if (F_dir.compare("b") == 0) {
			pos_a = rand()%n_box_a;
			pos_b = 0;
			pos_c = rand()%n_box_c;
		}
		
		else if (F_dir.compare("c") == 0) {
			pos_a = rand()%n_box_a;
			pos_b = rand()%n_box_b;
			pos_c = 0;
		}	
		
		mol = rand()%n_mol;
		
		//k_sum = 0.0;
		//
		//for(unsigned int jj=0; jj<neigh_label[frame][mol].size(); jj++) {
		//	k_sum = k_sum + k_inv[frame][mol][jj];
		//}
		
		// Find the number of the corresponding mini-grid
		for (unsigned int x=0; x<box_neigh_label.size(); x++){
			if (pos_a == box_a[x] && pos_b == box_b[x] && pos_c == box_c[x]){
				box = x;
				break;
			}
		}
	}
	while(neigh_label[frame][mol].size() == 0 || grid_occ[box][mol] == true);
	
	pos[0] = mol;
	pos[1] = box;
}

void Dispatch_Mol_begin_layer(int frame, vector< vector<bool> > grid_occ, int *pos){ 
	
	int pos_a=0, pos_b=0, pos_c=0, mol=0, box=0; 
	
	//double k_sum;
	
	do {
		if (F_dir.compare("a") == 0) {
			pos_a = 0;
			pos_b = rand()%n_box_b;
			pos_c = rand()%n_box_c;
		}
		
		else if (F_dir.compare("b") == 0) {
			pos_a = rand()%n_box_a;
			pos_b = 0;
			pos_c = rand()%n_box_c;
		}
		
		else if (F_dir.compare("c") == 0) {
			pos_a = rand()%n_box_a;
			pos_b = rand()%n_box_b;
			pos_c = 0;
		}	
		
		//mol = rand()%n_mol;
		mol = mol_layer[frame][layer][rand()%mol_layer[frame][layer].size()];
		
		//k_sum = 0.0;
		//
		//for(unsigned int jj=0; jj<neigh_label[frame][mol].size(); jj++) {
		//	k_sum = k_sum + k_inv[frame][mol][jj];
		//}
		
		// Find the number of the corresponding mini-grid
		for (unsigned int x=0; x<box_neigh_label.size(); x++){
			if (pos_a == box_a[x] && pos_b == box_b[x] && pos_c == box_c[x]){
				box = x;
				break;
			}
		}
	}
	while(neigh_label[frame][mol].size() == 0 || grid_occ[box][mol] == true);
	
	pos[0] = mol;
	pos[1] = box;
}
