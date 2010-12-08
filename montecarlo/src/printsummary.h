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

#ifndef _PRINTSUMMARY_H
#define _PRINTSUMMARY_H 1

// Print all the informations used to calculate k (before running any simulation)
void Print_Summary_Beginning(string output_folder);

// Print some information for each try of the MC simulation
void Print_Summary_Try(string output_folder, int i, int charge_try, double total_dist_try,\
																					double total_time_try);
																					
// Print some information at the end of each frame calculation
void Print_Summary_Frame(string output_folder, int i, double total_dist_try, double total_time_try,\
																					vector <double> mu_frame);
																					
// Print some information at the end of the simulation
void Print_Summary_Final(string output_folder, double mu_moy);
	
#endif // printsummary.h
