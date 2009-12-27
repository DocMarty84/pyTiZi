#!/usr/bin/env python

import sys
import numpy as np
import time
import copy
import os, sys, shutil

sys.path.append("lib")
import read_input_MD, molecular_system, list_manipulation, write_cluster_files, coord_conversion
from data_input import *

if __name__ == '__main__':
	""" 
		This program takes as an input the coordinates provided by CRYSTAL
		and reproduces them along several axis, in order 
		to get 2 entire	molecules.
	"""
	
	# First part: reproduces the coordinates in several directions, and
	# writes the file caca.arc you should open with VMD. With VMD, get the
	# index of each atom (do not forget to add 1 to the index).
	
	try:
		finput = open("coordinates.xyz", 'r')
	except:
		print "Could not open %s" % ("coord.xyz")
		sys. exit(1)
	
	line = finput.readline()
	words = line.split()
	n_frame = int(words[0])
	n_mol = int(words[1])
	n_atom = int(words[2])
	
	a = np.zeros((1), float)
	b = np.zeros((1), float)
	c = np.zeros((1), float)
	alpha = np.zeros((1), float)
	beta = np.zeros((1), float)
	gamma = np.zeros((1), float)
	
	line = finput.readline()
	words = line.split()
	a[0] = float(words[0])
	b[0] = float(words[1])
	c[0] = float(words[2])
	alpha[0] = float(words[3])
	beta[0] = float(words[4])
	gamma[0] = float(words[5])
	
	coord = molecular_system.MolecularSystem(n_frame, n_mol, [n_atom])
	
	for i in xrange(coord.n_atom[0]):
		line = finput.readline()
		words = line.split()
		
		coord.symbol[0][0][i] = words[2]
		coord.x[0,0,i] = float(words[3])
		coord.y[0,0,i] = float(words[4])
		coord.z[0,0,i] = float(words[5])
	
	finput.close()

	# Get the cell parameters as well as the number frame, molecules and atoms
	box = molecular_system.SimulationBox(n_frame, a, b, c, alpha, beta, gamma)
	box.Parameters_For_Orthogonalization()
	
	k = 0
	tmp = "%d ANTH\n" % (coord.n_atom[0])

	for i in [-1, 0, 1]:				
		for ii in [0, 1]:
			for iii in [0, 1]:
				for t in xrange(coord.n_atom[0]):
					k += 1
					Frac_Coord = coord_conversion.Cartesian_To_Fractional(coord.x[0,0,t], coord.y[0,0,t], coord.z[0,0,t], box)
					Frac_Coord[0] += i
					Frac_Coord[1] += ii
					Frac_Coord[2] += iii
					Cart_Coord = coord_conversion.Fractional_To_Cartesian(Frac_Coord[0], Frac_Coord[1], Frac_Coord[2], box)
					tmp += "%4d %3s %15.10f %15.10f %15.10f %6d\n" % (k, coord.symbol[0][0][t], Cart_Coord[0], Cart_Coord[1], Cart_Coord[2], t+1)
			
	try:
		finput = open("caca.arc", 'w')
	except:
		print "Could not open %s" % ("caca.arc")
		sys. exit(1)
		
	finput.write(tmp)
	finput.close
	
	# At this point, you should have the index of the atoms of both molecules.
	# Write these index in the 2 following list, and the program will print
	# only the good atoms.
	
	"""
	mol1=[196, 316, 312, 308, 302, 346, 342, 338, 266, 262, 258, 256, 204, 200, 464, 336, 332, 326, 370, 126, 286, 282, 280, 228]
	mol2=[241, 265, 261, 257, 255, 251, 247, 243, 267, 263, 259, 253, 249, 245, 269, 285, 281, 279, 275, 271, 287, 283, 277, 273]

	try:
		finput = open("caca.arc", 'r')
	except:
		print "Could not open %s" % ("caca.arc")
		sys. exit(1)
	
	tmp = 'mol1\n'
	for line in finput:
		words = line.split()
		
		if int(words[0]) in mol1:
			tmp += line

	finput.close

	try:
		finput = open("caca.arc", 'r')
	except:
		print "Could not open %s" % ("caca.arc")
		sys. exit(1)
		
	tmp += 'mol2\n'
	for line in finput:
		words = line.split()
		
		if int(words[0]) in mol2:
			tmp += line
			
	finput.close

	print tmp
	"""

	sys.exit(0)
