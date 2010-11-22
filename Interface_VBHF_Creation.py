#!/usr/bin/env python

import sys
import numpy as np
import time
import copy
import os, sys, shutil

sys.path.append("lib")
import read_input_MD, molecular_system, list_manipulation, write_cluster_files, coord_conversion
from data_input import *

TEMPERATURE = 300						# Temperature in K
H_BAR = 3.990312682e+13/(2*np.pi)		# Planck constant in atomic mass units and angstrom^2
H_BAR_EV = 6.58211899e-16				# Planck constant in eV
K_B = 8.314472477e+23					# Boltmann constant in atomic mass units and angstrom^2
CM1_TO_HZ = 2.99792458e+10*(2*np.pi)	# conversion between cm-1 and Hz, multiplied by 2pi for pulsation


if __name__ == '__main__':
	""" 
		This program generates the quantum sampling structures from the
		eigenvectors of the Hessian matrix.
		At the moment, it only accepts mass-weighted eigenvectors.
	"""
	t1 = time.clock()
	
	filename_base = "vertical300K"
	file_pdb = "%s.pdb" % (filename_base)
	
	# Get the cell parameters as well as the number frame, molecules and atoms
	(n_frame, n_mol, n_atom, a, b, c, alpha, beta, gamma) = read_input_MD.Read_PDB_File_First(file_pdb)
	box = molecular_system.SimulationBox(n_frame, a, b, c, alpha, beta, gamma)
	box.Parameters_For_Orthogonalization()
	
	#print n_frame, n_mol, n_atom
	
	# Get the coordinates of the system
	coord = molecular_system.MolecularSystem(n_frame, n_mol, n_atom)
	read_input_MD.Read_PDB_File_Second(file_pdb, coord)
	#print coord
	
	write_cluster_files.CreateVBHFInput_Interface(coord, box)

	t2 = time.clock()
	print t2-t1
	
	sys.exit(0)
