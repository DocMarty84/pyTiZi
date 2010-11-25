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

#ifndef _CONSTANTS_H
#define _CONSTANTS_H 1

#ifdef INSIDE_MC_BKL
#define EXTERN
#else
#define EXTERN extern
#endif

// Physical contants in the correct units
EXTERN const double K_BOLTZ		= 8.617343e-5;			// Boltzmann const (eV)
EXTERN const double H_BAR		= 6.58211899e-16;		// h/2PI in eV
EXTERN const double EPSILON_0	= (0.5)*(137.035999084)* (1.0/4.13566733e-15)\
								  *(1.0/299792458e10); 	// Vacuum permittivity 
														// e/(V.Ang)
EXTERN const double PI 		= 3.14159265358979323846;	// Pi

#endif
