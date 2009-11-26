#!/usr/bin/env python

import sys
import numpy
import time
import copy

sys.path.append("lib")
import read_input_MD, molecular_system, list_manipulation, write_cluster_files
from data_input import *


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
	
	(n_frame, n_mol, n_atom, a, b, c, alpha, beta, gamma) = read_input_MD.Read_TINKER_add_File(file_add)
	box = molecular_system.SimulationBox(n_frame, a, b, c, alpha, beta, gamma)
	
	qs_coord = molecular_system.MolecularSystem(n_frame, n_mol, [n_atom])
	read_input_MD.Read_TINKER_arc_File(file_xyz, qs_coord)
#	print qs_coord
	
	qs_eigen = molecular_system.Eigenvector_Matrix(qs_coord.n_atom[0]*3)
	read_input_MD.Read_Tinker_vec_File(file_vec, qs_coord, qs_eigen)
#	print qs_eigen
#	print sum(numpy.ravel(qs_eigen.vec_matrix[:,:])*numpy.ravel(qs_eigen.vec_matrix[:,:]))
	
	qs_eigen.getI()
#	print qs_eigen.vec_matrix*qs_eigen.vec_matrix_I

	qs_eigen.getRefCoord(qs_coord)
#	print qs_eigen.cart_coord_matrix

	qs_eigen.getNormalCoord()
#	print qs_eigen.normal_coord_matrix

	tmp = ''
	Q_start = qs_eigen.normal_coord_matrix
	for d in xrange(1,10000,1000):
		Q = copy.copy(Q_start)
		X = qs_eigen.vec_matrix_I*Q
		
		j = 0
		X2 = numpy.zeros((qs_coord.n_mol*qs_coord.n_atom), float)
		mu_up = 0.0
		mu_down = 0.0
		for i in xrange(qs_coord.n_mol*qs_coord.n_atom):
			a = i%qs_coord.n_mol
			b = i%qs_coord.n_atom[i%qs_coord.n_mol]
					
			X2[i] = X[j, 0]*X[j, 0] + X[j+1, 0]*X[j+1, 0] + X[j+2, 0]*X[j+2, 0]
			mu_up = mu_up + qs_coord.atomic_mass[0, a, b] * X2[i]
			mu_down = qs_coord.atomic_mass[0, a, b] * X2[i]

			j += 3

		mu = numpy.sqrt(mu_up) / numpy.sqrt(mu_down)
#		print mu
			
		Q[0,0] = Q[0,0]*d*numpy.sqrt(2/mu)
#		print Q
		X = qs_eigen.vec_matrix_I*Q
#		print X
		
		tmp += '%5d %s\n' % (qs_coord.n_mol*qs_coord.n_atom, filename_base)
		j = 0
		for i in xrange(qs_coord.n_mol*qs_coord.n_atom):
			a = i%qs_coord.n_mol
			b = i%qs_coord.n_atom[i%qs_coord.n_mol]
			tmp += '%5d %2s %12.6f %12.6f %12.6f 0\n' % (i + 1, qs_coord.symbol[0][a][b], X[j, 0], X[j+1, 0], X[j+2, 0])
			j += 3
		
	print tmp

	t2 = time.clock()
	print t2-t1
