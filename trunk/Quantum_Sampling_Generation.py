#!/usr/bin/env python

import sys
import numpy as np
import time
import copy
import os, sys

sys.path.append("lib")
import read_input_MD, molecular_system, list_manipulation, write_cluster_files, coord_conversion
from data_input import *

TEMPERATURE = 300						# Temperature in K
H_BAR = 3.990312682e+13/(2*np.pi)		# Planck constant in atomic mass units and angstrom^2
H_BAR_EV = 6.58211899e-16				# Planck constant in eV
K_B = 8.314472477e+23					# Boltmann constant in atomic mass units and angstrom^2
CM1_TO_HZ = 2.99792458e+10*(2*np.pi)	# conversion between cm-1 and Hz, multiplied by 2pi for pulsation

DELTA_E_MAX = 0.1*K_B*TEMPERATURE

CHECK = False
MODE = 40

if __name__ == '__main__':
	""" 
		This program generates the quantum sampling structures from the
		eigenvectors of the Hessian matrix.
		At the moment, it only accepts mass-weighted eigenvectors.
	"""
	t1 = time.clock()
	
	filename_base = "anthracene_UNIT_CELL_MIN"
	file_vec = "%s.vec" % (filename_base)
	file_xyz = "%s.xyz" % (filename_base)
	file_add = "%s.add" % (filename_base)
	
	# Get the cell parameters as well as the number frame, molecules and atoms
	(n_frame, n_mol, n_atom, a, b, c, alpha, beta, gamma) = read_input_MD.Read_TINKER_add_File(file_add)
	box = molecular_system.SimulationBox(n_frame, a, b, c, alpha, beta, gamma)
	box.Parameters_For_Orthogonalization()
	
	# Get the coordinates of the system
	qs_coord = molecular_system.MolecularSystem(n_frame, n_mol, [n_atom])
	read_input_MD.Read_TINKER_arc_File(file_xyz, qs_coord)
#	print qs_coord
	
	# Get the eigenvectors of the system. After this part, the eigenvectors
	# must be non mass-weighted and normalized
	qs_eigen = molecular_system.Eigenvector_Matrix(qs_coord.n_mol*qs_coord.n_atom[0]*3)
	read_input_MD.Read_TINKER_vec_File(file_vec, qs_coord, qs_eigen)
#	print qs_eigen
#	print sum(np.ravel(qs_eigen.vec_matrix[:,:])*np.ravel(qs_eigen.vec_matrix[:,:]))

	# Copy the coordinates of the equilibrium system
	qs_eigen.getRefCoord(qs_coord)
#	print qs_eigen.cart_coord_matrix

	# Get the mass of each atoms
	qs_eigen.getMasses(qs_coord)
#	print qs_eigen.mass_matrix

	# Converts the eigenvectors read in the output file in non mass-weighted
	# and normalized eigenvectors
	qs_eigen.ConvertEigenVec("tinker")
#	print sum(np.ravel(qs_eigen.vec_matrix[1,:])*np.ravel(qs_eigen.vec_matrix[1,:]))

	# Calculates Normal Coordinates, on the basis of the 
	# mass-weighted eigenvectors
#	qs_eigen.getNormalCoord("tinker")
#	print qs_eigen.normal_coord_matrix

	# Calculates the inverse matrix of both eigenvectors matrix
#	qs_eigen.getI()
#	print qs_eigen.vec_matrix*qs_eigen.vec_matrix_I
#	print qs_eigen.vec_matrix_MW*qs_eigen.vec_matrix_MW_I

	X_ref = copy.copy(qs_eigen.cart_coord_matrix)
	T_ref = copy.copy(qs_eigen.vec_matrix)
		
	for mode in xrange(0, qs_eigen.n_modes, 1):
		if qs_eigen.freq[mode] > 0.0:
#	for mode in xrange(12, 13, 1):
			# ===============================
			# Calculation of the reduced mass
			# ===============================
			j = 0
			mu_up = 0.0
			mu_down = 0.0
			for i in xrange(qs_coord.n_mol*qs_coord.n_atom[0]):
				a = i%qs_coord.n_mol
				b = i%qs_coord.n_atom[i%qs_coord.n_mol]

#				t2 = 0.0
				for k in xrange(3):
					t2 = np.power(qs_eigen.vec_matrix[mode,j+k], 2)
					mu_down = mu_down + qs_coord.atomic_mass[0, a, b] * t2
					mu_up = mu_up + t2

				j += 3

			mu_I = mu_up / mu_down
			mu = 1/mu_I
			
			d_max = np.sqrt(2*DELTA_E_MAX/(H_BAR*CM1_TO_HZ*qs_eigen.freq[mode]))
			
			if MODE==mode+1 and CHECK:
				print H_BAR_EV*CM1_TO_HZ*qs_eigen.freq[mode]/2
			
			tmp = ''
			p = False
			
			for d_step in xrange(-10,11,2):
				d = (d_max)*(d_step/10.0)
				
				if MODE==mode+1 and CHECK:
					print d
				
				X = copy.copy(X_ref)
				T = copy.copy(T_ref)
#				print sum(np.ravel(T[:,:])*np.ravel(T[:,:]))
				
				# ==========================================
				# Modify the eigenvectors of the normal mode
				# ==========================================
				
				# Apply the displacement
				T[mode,:] = np.sqrt(H_BAR/(mu*CM1_TO_HZ*qs_eigen.freq[mode]))*T[mode,:]*d

				X_T = X.getT()
				X_T = X_T + T[mode,:]
				X = X_T.getT()
				
				# Create Tinker file for visualization of the mode
				if d == d_max:
					p = True
				tmp = write_cluster_files.CreateTINKERVisualization(X, qs_coord, mode, tmp, filename_base, p)
				
				# Create Normal Modes files like Sigi
#				write_cluster_files.CreateNormalMode(X, qs_coord, mode, d)
				
				# Create VBHF input files
#				write_cluster_files.CreateVBHFInput(X, qs_coord, box, mode, d)

	# Create a script which will create all the needed pbs
#	write_cluster_files.ScriptVBHFLaunch("/home/nmartine/VBHF/QS_Tinker")
	
	t2 = time.clock()
	print t2-t1
	
	sys.exit(0)
