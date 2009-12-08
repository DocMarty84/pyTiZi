#!/usr/bin/env python

import sys
import numpy as np
import time
import copy
import os, sys

sys.path.append("lib")
import read_input_MD, molecular_system, list_manipulation, write_cluster_files, coord_conversion
from data_input import *

h_planck = 3.990312682e+13	# Planck constant in atomic mass units and angstrom^2
cm1_TO_Hz = 2.99792458e+10	# conversion between cm-1 and Hz

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
#	print sum(np.ravel(qs_eigen.vec_matrix[:,:])*np.ravel(qs_eigen.vec_matrix[:,:]))

	# Calculates Normal Coordinates, on the basis of the 
	# mass-weighted eigenvectors
#	qs_eigen.getNormalCoord("tinker")
#	print qs_eigen.normal_coord_matrix

	# Calculates the inverse matrix of both eigenvectors matrix
#	qs_eigen.getI()
#	print qs_eigen.vec_matrix*qs_eigen.vec_matrix_I
#	print qs_eigen.vec_matrix_MW*qs_eigen.vec_matrix_MW_I

	X_ref   = copy.copy(qs_eigen.cart_coord_matrix)
	T_ref   = copy.copy(qs_eigen.vec_matrix)
	
	try:
		os.mkdir("normal_modes")
	except:
		pass
		
#	for mode in xrange(4, qs_eigen.n_modes, 1):
#		print qs_eigen.freq[mode]
	for mode in xrange(12, 13, 1):
		# =================================
		# Calculation of the "reduced mass"
		# =================================
		j = 0
		mu_up = 0.0
		mu_down = 0.0
		for i in xrange(qs_coord.n_mol*qs_coord.n_atom[0]):
			a = i%qs_coord.n_mol
			b = i%qs_coord.n_atom[i%qs_coord.n_mol]

			t2 = 0.0
			for k in xrange(3):
				t2 += np.power(qs_eigen.vec_matrix[mode,j+k], 2)
			mu_up = mu_up + qs_coord.atomic_mass[0, a, b] * t2
			mu_down = mu_down + t2

			j += 3

		mu = mu_up / mu_down
	
		dir = "normal_modes/normal_mode_%d" % (mode+1)
		try:
			os.mkdir(dir)
		except:
			pass
			
		for d in xrange(-5,6,1):
			tmp = ''
			d = d/10.0
			X = copy.copy(X_ref)
			T = copy.copy(T_ref)
#			print sum(np.ravel(T[:,:])*np.ravel(T[:,:]))
			
			# ==========================================
			# Modify the eigenvectors of the normal mode
			# ==========================================
			
			# Apply the displacement
			T[mode,:] = np.sqrt(h_planck/(mu*cm1_TO_Hz*qs_eigen.freq[mode]))*T[mode,:]*d

			X_T = X.getT()
			X_T = X_T + T[mode,:]
			X = X_T.getT()
			
#			tmp += '%5d %s\n' % (qs_coord.n_mol*qs_coord.n_atom[0], filename_base)
#			j = 0
#			for i in xrange(qs_coord.n_mol*qs_coord.n_atom[0]):
#				a = i%qs_coord.n_mol
#				b = i%qs_coord.n_atom[i%qs_coord.n_mol]
#				tmp += '%5d %2s %12.6f %12.6f %12.6f 0\n' % (i + 1, qs_coord.symbol[0][a][b], X[j, 0], X[j+1, 0], X[j+2, 0])
#				j += 3
				
			j = 0
			for i in xrange(qs_coord.n_mol*qs_coord.n_atom[0]):
				a = i%qs_coord.n_mol
				b = i%qs_coord.n_atom[i%qs_coord.n_mol]
				tmp += '%2s %12.6f %12.6f %12.6f\n' % (qs_coord.symbol[0][a][b], X[j, 0], X[j+1, 0], X[j+2, 0])
				j += 3
			
			if d < 0.0:
				name = "%s/result-%d-minus-%.1f.dat" % (dir, mode+1, abs(d))
			else:
				name = "%s/result-%d-plus-%.1f.dat" % (dir, mode+1, d)
				
			try:
				foutput = open(name, 'w')
			except:
				print "Could not open %s" % (name)
				sys.exit(1)
				
			foutput.write(tmp)
			foutput.close()	
			
######################################

			for charge in [-1, 0, 1]:
				try:
					dir_all="VBHF/all_%d" % (charge)
					os.makedirs(dir_all)
				except:
					pass
				try:
					dir_mono="VBHF/mono_%d" % (charge)
					os.makedirs(dir_mono)
				except:
					pass
							
				print "Calculating normal mode %d of charge %d" %(mode+1, charge)
					
				# All cluster
				if d < 0.0:
					name = "%s/result-%d-minus-%.1f.dat" % (dir_all, mode+1, abs(d))
				else:
					name = "%s/result-%d-plus-%.1f.dat" % (dir_all, mode+1, d)
					
				foutput = open(name, 'w')
				
				if foutput:
					tmp = ''
					tmp = "AM1 1SCF VBHF\n\n"
					tmp += "Xx        0.0000 1     0.0000 1     0.0000 1\n"
					tmp += "Xx        0.0000 1     0.0000 1     0.0000 1\n"
					tmp += "Xx        1.0000 1     0.0000 1     0.0000 1\n"
					tmp += "Xx        0.0000 1     1.0000 1     0.0000 1\n"
					tmp += "Xx        0.0000 1     0.0000 1     1.0000 1\n" 

					k = 1
					for a in [0, -1, 1]:
						for b in [0, -1, 1]:
							j = 0
							for i in xrange(qs_coord.n_atom[0]):
								Frac_Coord = coord_conversion.Cartesian_To_Fractional(X[j, 0], X[j+1, 0], X[j+2, 0], box)
								Frac_Coord[0] += a
								Frac_Coord[1] += b
								Cart_Coord = coord_conversion.Fractional_To_Cartesian(Frac_Coord[0], Frac_Coord[1], Frac_Coord[2], box)
								tmp += "%4s %12f 1 %12f 1 %12f 1\n" % (qs_coord.symbol[0][0][i], Cart_Coord[0], Cart_Coord[1], Cart_Coord[2])
								#tmp += "%5d %4s %12f %12f %12f 0\n" % (k, qs_coord.symbol[0][1][i], Cart_Coord[0], Cart_Coord[1], Cart_Coord[2])
								j += 3
								k += 1
							for i in xrange(qs_coord.n_atom[1]):
								Frac_Coord = coord_conversion.Cartesian_To_Fractional(X[j, 0], X[j+1, 0], X[j+2, 0], box)
								Frac_Coord[0] += a
								Frac_Coord[1] += b
								Cart_Coord = coord_conversion.Fractional_To_Cartesian(Frac_Coord[0], Frac_Coord[1], Frac_Coord[2], box)
								tmp += "%4s %12f 1 %12f 1 %12f 1\n" % (qs_coord.symbol[0][1][i], Cart_Coord[0], Cart_Coord[1], Cart_Coord[2])
								#tmp += "%5d %4s %12f %12f %12f 0\n" % (k, qs_coord.symbol[0][1][i], Cart_Coord[0], Cart_Coord[1], Cart_Coord[2])
								j += 3
								k += 1
								
					for a in [-2, -1, 0, 1]:
						j = 0
						for i in xrange(qs_coord.n_atom[1]):
							Frac_Coord = coord_conversion.Cartesian_To_Fractional(X[(qs_coord.n_atom[0]*3)+j, 0], X[(qs_coord.n_atom[0]*3)+j+1, 0], X[(qs_coord.n_atom[0]*3)+j+2, 0], box)
							Frac_Coord[0] += a
							Frac_Coord[1] += -2
							Cart_Coord = coord_conversion.Fractional_To_Cartesian(Frac_Coord[0], Frac_Coord[1], Frac_Coord[2], box)
							tmp += "%4s %12f 1 %12f 1 %12f 1\n" % (qs_coord.symbol[0][0][i], Cart_Coord[0], Cart_Coord[1], Cart_Coord[2])
							#tmp += "%5d %4s %12f %12f %12f 0\n" % (k, qs_coord.symbol[0][0][i], Cart_Coord[0], Cart_Coord[1], Cart_Coord[2])
							j += 3
							k += 1
								
					for b in [-1, 0, 1]:
						j = 0
						for i in xrange(qs_coord.n_atom[1]):
							Frac_Coord = coord_conversion.Cartesian_To_Fractional(X[(qs_coord.n_atom[0]*3)+j, 0], X[(qs_coord.n_atom[0]*3)+j+1, 0], X[(qs_coord.n_atom[0]*3)+j+2, 0], box)
							Frac_Coord[0] += -2
							Frac_Coord[1] += b
							Cart_Coord = coord_conversion.Fractional_To_Cartesian(Frac_Coord[0], Frac_Coord[1], Frac_Coord[2], box)
							tmp += "%4s %12f 1 %12f 1 %12f 1\n" % (qs_coord.symbol[0][0][i], Cart_Coord[0], Cart_Coord[1], Cart_Coord[2])
							#tmp += "%5d %4s %12f %12f %12f 0\n" % (k, qs_coord.symbol[0][0][i], Cart_Coord[0], Cart_Coord[1], Cart_Coord[2])
							j += 3
							k += 1
				
					tmp += "$$VBHF\n"
					tmp += "%d %d AM1 OMF-OPT\n" % (qs_coord.n_atom[0], charge)
					for x in xrange(qs_coord.n_atom[0]-1):
						tmp += "%d 0 AM1 OMF-OPT\n"	 % qs_coord.n_atom[0]	

					foutput.write(tmp)		
					foutput.close() 
		
		
				# Molecule alone
				if d < 0.0:
					name = "%s/result-%d-minus-%.1f.dat" % (dir_mono, mode+1, abs(d))
				else:
					name = "%s/result-%d-plus-%.1f.dat" % (dir_mono, mode+1, d)
					
				foutput = open(name, 'w')
				
				if foutput:
					tmp = ''
					tmp = "AM1 1SCF VBHF\n\n"
					tmp += "Xx        0.0000 1     0.0000 1     0.0000 1\n"
					tmp += "Xx        0.0000 1     0.0000 1     0.0000 1\n"
					tmp += "Xx        1.0000 1     0.0000 1     0.0000 1\n"
					tmp += "Xx        0.0000 1     1.0000 1     0.0000 1\n"
					tmp += "Xx        0.0000 1     0.0000 1     1.0000 1\n"

					for i in xrange(qs_coord.n_atom[0]):
						tmp += "%4s %12f 1 %12f 1 %12f 1\n" % (qs_coord.symbol[0][1][i], X[i, 0], X[i+1, 0], X[i+2, 0])
					
					tmp += "$$VBHF\n"
					tmp += "%d %d AM1 OMF-OPT\n" % (qs_coord.n_atom[0], charge)	

					foutput.write(tmp)		
					foutput.close()


######################################

	t2 = time.clock()
	print t2-t1
	
	sys.exit(0)
