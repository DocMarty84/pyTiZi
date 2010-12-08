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

// C Libraries
#include <stdlib.h>

// Local libraries
#include "variables.h"
#include "coordinates.h"

using namespace std;

// ===========================================================================================================
// ---------------------------------------------- Grid building ----------------------------------------------
// ===========================================================================================================

// Find the layer for each molecule. Need to be modified.
void Find_Layer(bool print_results) {
	min_layer[0]=-100.0; max_layer[0]=-40.0;
	min_layer[1]=-35.0; max_layer[1]=-25.0;
	min_layer[2]=-22.0; max_layer[2]=-12.0;
	min_layer[3]=-10.0; max_layer[3]=10.0;
	
	for (int i=0; i<n_frame; i++){
		list_layer.push_back(vector< int > ());
		mol_layer.push_back(vector< vector<int> > ());	
		for (int k=0; k<n_layer; k++){	
			mol_layer[i].push_back(vector< int > ());
		}
	}

	for (int i=0; i<n_frame; i++){
		for (int ii=0; ii<n_mol; ii++){
			for (int k=0; k<n_layer; k++){
				if(CM_z[i][ii] < max_layer[k] && CM_z[i][ii] > min_layer[k]){
					mol_layer[i][k].push_back(ii);
					list_layer[i].push_back(k);
				}
			}
		}
	}
	
	if(print_results){
		for (int i=0; i<n_frame; i++){
			cout << endl;
			cout << "Frame " << i+1 << endl;
			for (int k=0; k<n_layer; k++){
				cout << endl;
				cout << "Layer " << k << ", with " << mol_layer[i][k].size() <<	" molecules" << endl;
				for (unsigned int ii=0; ii<mol_layer[i][k].size(); ii++){
					cout << mol_layer[i][k][ii] << " ";
				}
			}
			cout << endl;
		}
	}
}

// Build the grid of molecules, and define the PBC for the whole system
void Build_Grid(bool print_results) {
	
	int n_grid = n_box_a*n_box_b*n_box_c;
	
	// Set up the "coordinates" of each box
	for (int a=0; a<n_box_a; a++){
		for (int b=0; b<n_box_b; b++){
			for (int c=0; c<n_box_c; c++){
				box_a.push_back(a);
				box_b.push_back(b);
				box_c.push_back(c);

			}
		}
	}
	
	// Set up the neighbors of each box
	for (int x=0; x<n_grid; x++){
		box_neigh_label.push_back( vector< int > () );
		box_neigh_a.push_back( vector< int > () );
		box_neigh_b.push_back( vector< int > () );
		box_neigh_c.push_back( vector< int > () );
		
		for (int y=0; y<n_grid; y++){
			// Finite size system
			if ((abs(box_a[y]-box_a[x]) < 2 && abs(box_b[y]-box_b[x]) < 2 && abs(box_c[y]-box_c[x]) < 2) &&\
				(abs(box_a[y]-box_a[x]) != 0 || abs(box_b[y]-box_b[x]) != 0 || abs(box_c[y]-box_c[x]) != 0)){
					
				box_neigh_a[x].push_back( box_a[y]-box_a[x] );
				box_neigh_b[x].push_back( box_b[y]-box_b[x] );
				box_neigh_c[x].push_back( box_c[y]-box_c[x] );
				box_neigh_label[x].push_back( y );
			}

			// Create an infinite system in the case of an anisotropy calculation
			if (anisotropy && box_a[x] == 0 && box_a[y] == n_box_a-1 && abs(box_b[y]-box_b[x]) < 2 &&\
																				abs(box_c[y]-box_c[x]) < 2) {
						
				box_neigh_a[x].push_back( -1 );
				box_neigh_b[x].push_back( box_b[y]-box_b[x] );
				box_neigh_c[x].push_back( box_c[y]-box_c[x] );
				box_neigh_label[x].push_back( y );
			}

			if (anisotropy && box_b[x] == 0 && box_b[y] == n_box_b-1 && abs(box_a[y]-box_a[x]) < 2 &&\
																				abs(box_c[y]-box_c[x]) < 2) {
						
				box_neigh_b[x].push_back( -1 );
				box_neigh_a[x].push_back( box_a[y]-box_a[x] );
				box_neigh_c[x].push_back( box_c[y]-box_c[x] );
				box_neigh_label[x].push_back( y );
			}

			if (anisotropy && box_c[x] == 0 && box_c[y] == n_box_c-1 && abs(box_a[y]-box_a[x]) < 2 &&\
																				abs(box_b[y]-box_b[x]) < 2) {
						
				box_neigh_c[x].push_back( -1 );
				box_neigh_a[x].push_back( box_a[y]-box_a[x] );
				box_neigh_b[x].push_back( box_b[y]-box_b[x] );
				box_neigh_label[x].push_back( y );
			}

			if (anisotropy && box_a[x] == n_box_a-1 && box_a[y] == 0 && abs(box_b[y]-box_b[x]) < 2 &&\
																				abs(box_c[y]-box_c[x]) < 2) {
				box_neigh_a[x].push_back( 1 );
				box_neigh_b[x].push_back( box_b[y]-box_b[x] );
				box_neigh_c[x].push_back( box_c[y]-box_c[x] );
				box_neigh_label[x].push_back( y );
			}

			if (anisotropy && box_b[x] == n_box_b-1 && box_b[y] == 0 && abs(box_a[y]-box_a[x]) < 2 &&\
																				abs(box_c[y]-box_c[x]) < 2) {
				box_neigh_b[x].push_back( 1 );
				box_neigh_a[x].push_back( box_a[y]-box_a[x] );
				box_neigh_c[x].push_back( box_c[y]-box_c[x] );
				box_neigh_label[x].push_back( y );
			}

			if (anisotropy && box_c[x] == n_box_c-1 && box_c[y] == 0 && abs(box_a[y]-box_a[x]) < 2 &&\
																				abs(box_b[y]-box_b[x]) < 2) {
				box_neigh_c[x].push_back( 1 );
				box_neigh_a[x].push_back( box_a[y]-box_a[x] );
				box_neigh_b[x].push_back( box_b[y]-box_b[x] );
				box_neigh_label[x].push_back( y );
			}
		}
	}
	
	// Set up the grid properties: coordinates, probability and energy of each site
	double *CM_Cart, *CM_Frac;
	CM_Cart = new double[3];
	CM_Frac = new double[3];
	
	for (int i=0; i<n_frame; i++){
		grid_probability.push_back( vector< vector< vector<double> > > () );
		grid_x.push_back( vector< vector<double> > () );
		grid_y.push_back( vector< vector<double> > () ); 
		grid_z.push_back( vector< vector<double> > () );
		grid_E_0.push_back( vector< vector<double> > () );
		grid_E_1.push_back( vector< vector<double> > () );
		
		for (int x=0; x<n_grid; x++){
			
			grid_probability[i].push_back( vector< vector<double> > () );
			grid_x[i].push_back( vector<double> () );
			grid_y[i].push_back( vector<double> () ); 
			grid_z[i].push_back( vector<double> () );
			grid_E_0[i].push_back( vector<double> () );
			grid_E_1[i].push_back( vector<double> () );
			
			for (int ii=0; ii<n_mol; ii++){
				
				CM_Cart[0] = CM_x[i][ii];
				CM_Cart[1] = CM_y[i][ii];
				CM_Cart[2] = CM_z[i][ii];
				Cartesian_To_Fractional(CM_Cart, CM_Frac, i);
				
				CM_Frac[0] += double(box_a[x]);
				CM_Frac[1] += double(box_b[x]);
				CM_Frac[2] += double(box_c[x]);
				Fractional_To_Cartesian(CM_Frac, CM_Cart, i);
		
				grid_probability[i][x].push_back( vector<double> () );
				grid_x[i][x].push_back( CM_Cart[0] );
				grid_y[i][x].push_back( CM_Cart[1] ); 
				grid_z[i][x].push_back( CM_Cart[2] );
				grid_E_0[i][x].push_back( E_0[0][ii] );
				grid_E_1[i][x].push_back( E_1[0][ii] );
				
				// Set up a probability for each charge
				for (unsigned int charge_i = 0; charge_i < n_charges; charge_i++){
					grid_probability[i][x][ii].push_back( 0.0 );
				}
			}
		}
	}
	
	delete [] CM_Cart;
	delete [] CM_Frac;
	
	if (print_results){
		
		for (unsigned int x=0; x<box_neigh_label.size(); x++){
			cout << "mini-grid " << x << endl;			
			for (unsigned int xx=0; xx<box_neigh_label[x].size(); xx++){
				cout << box_neigh_label[x][xx] << " " << box_neigh_a[x][xx] << " " << box_neigh_b[x][xx] <<\
																			" " << box_neigh_c[x][xx] << endl;
			}
		}
		
		cout << "-----------------------------------------------------" << endl;
		
		for (int i=0; i<n_frame; i++){
			cout << endl << "Frame " << i << endl;
			cout << "Box Mol X Y Z E_0 E_1" << endl;
			
			for (int x=0; x<n_grid; x++){
				
				for (int ii=0; ii<n_mol; ii++){
					cout << x << " " << ii << " " << grid_x[i][x][ii] << " " <<	grid_y[i][x][ii] << " " <<\
						grid_z[i][x][ii] << " " << grid_E_0[i][x][ii] << " " << grid_E_1[i][x][ii] << endl;
				}
			}
		}
	}
}
