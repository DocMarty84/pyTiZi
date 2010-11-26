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

#ifndef _READ_H
#define _READ_H 1

// Read .mc file
void Read_MC(string input_file, string input_folder, bool print_results);

// Read .cell file
void Read_CELL(string input_file, string input_folder, bool print_results);

// Read .cm part
void Read_CM(string input_file, string input_folder, bool print_results);

// Read e_av part (average energies)
void Read_E_av(string input_file, string input_folder, bool print_results);


#endif // read.h
