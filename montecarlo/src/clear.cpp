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
#include <limits>
#include <vector>

// Local libraries
#include "variables.h"

using namespace std;

// =============================================================================
// --------------------------------- Clear part --------------------------------
// =============================================================================

void Clear_All(){
	a.clear(); b.clear(); c.clear();
	alpha_deg.clear(); beta_deg.clear(); gamma_deg.clear();
	vol_box.clear(); 
	temp_alpha_cos.clear(); temp_beta_sin.clear(); temp_beta_cos.clear(); 
	temp_gamma_sin.clear(); temp_gamma_cos.clear(); temp_beta_term.clear(); 
	temp_gamma_term.clear(); 

	mol_label.clear();
	CM_x.clear(); CM_y.clear(); CM_z.clear(); 
	E_0.clear(); E_1.clear();
	
	chrg_E_electrostatic.clear(); chrg_E_0.clear(); chrg_E_1.clear();

	box_a.clear(); box_b.clear(); box_c.clear(); 
	box_neigh_a.clear(); box_neigh_b.clear(); box_neigh_c.clear();
	box_neigh_label.clear(); 
	grid_occ.clear();
	grid_probability.clear();
	grid_x.clear(); grid_y.clear(); grid_z.clear();
	grid_E_0.clear(); grid_E_1.clear();

	neigh_label.clear();
	d_x.clear(); d_y.clear(); d_z.clear();
	dE.clear();
	J_H.clear(); J_L.clear();
	neigh_jump_vec_a.clear(); 
	neigh_jump_vec_b.clear(); 
	neigh_jump_vec_c.clear(); 
	k.clear(); k_inv.clear();

	F_x_list.clear(); F_y_list.clear(); F_z_list.clear();
	F_angle_list.clear();

	MLJ_CST3.clear();

	min_layer.clear(); max_layer.clear();
	mol_layer.clear(); list_layer.clear();
}
