#!/usr/bin/env python

import string
import os, shutil, sys
import math
import numpy as np
import scipy as sp
import time
import copy

K_B = 8.617343e-5		# Boltzmann constant in eV

R_MIN = 0.90			# Minimum R^2 value to consider a mode as linear
N_TEMP = 500			# Temperature max for calculating sigma (starting at 1K)
FILE = "AM1"

#=====================================================================
#-------------------------Class Definition----------------------------
#=====================================================================

class VBHF_Data(object):
	""" Class for the VBHF data"""
	def __init__(self, n_modes=1, n_displ=1):
		self.n_modes = n_modes	# Number of vibrational modes
		self.n_displ = n_displ	# Number of displacements, including equilibrium geometry
		
		self.displ = np.zeros((n_modes, n_displ), float)			# Displacement in normal modes coordinates
		self.freq = np.zeros((n_modes), float)	
		self.energy = np.zeros((n_modes), float)	
		
		self.E_cluster_0 = np.zeros((n_modes, n_displ), float)	# VBHF energy for neutral cluster
		self.E_cluster_1 = np.zeros((n_modes, n_displ), float)	# VBHF energy for cluster with charge=1
		self.E_cluster_m1 = np.zeros((n_modes, n_displ), float)	# VBHF energy for cluster with charge=-1

		self.E_alone_0 = np.zeros((n_modes, n_displ), float)		# VBHF energy for alone neutral molecule
		self.E_alone_1 = np.zeros((n_modes, n_displ), float)		# VBHF energy for alone molecule with charge=1
		self.E_alone_m1 = np.zeros((n_modes, n_displ), float)	# VBHF energy for alone molecule with charge=-1

		self.IP_cluster = np.zeros((n_modes, n_displ), float)	# VBHF energy for alone molecule with charge=1
		self.EA_cluster = np.zeros((n_modes, n_displ), float)	# VBHF energy for alone molecule with charge=-1
		self.IP_alone = np.zeros((n_modes, n_displ), float)		# VBHF energy for alone molecule with charge=1
		self.EA_alone = np.zeros((n_modes, n_displ), float)		# VBHF energy for alone molecule with charge=-1
		self.P_plus = np.zeros((n_modes, n_displ), float)		# VBHF energy for alone molecule with charge=1
		self.P_minus = np.zeros((n_modes, n_displ), float)		# VBHF energy for alone molecule with charge=-1	

		self.IP_cluster_av = np.zeros((n_modes), float)	
		self.EA_cluster_av = np.zeros((n_modes), float)
		self.IP_alone_av = np.zeros((n_modes), float)
		self.EA_alone_av = np.zeros((n_modes), float)
		self.P_plus_av = np.zeros((n_modes), float)
		self.P_minus_av = np.zeros((n_modes), float)	
		
		self.IP_cluster_var = np.zeros((n_modes), float)	
		self.EA_cluster_var = np.zeros((n_modes), float)
		self.IP_alone_var = np.zeros((n_modes), float)
		self.EA_alone_var = np.zeros((n_modes), float)
		self.P_plus_var = np.zeros((n_modes), float)
		self.P_minus_var = np.zeros((n_modes), float)
	
		self.IP_cluster_dE = np.zeros((n_modes), float)	
		self.EA_cluster_dE = np.zeros((n_modes), float)
		self.IP_alone_dE = np.zeros((n_modes), float)
		self.EA_alone_dE = np.zeros((n_modes), float)
		self.P_plus_dE = np.zeros((n_modes), float)
		self.P_minus_dE = np.zeros((n_modes), float)	
		
		self.IP_cluster_R = np.zeros((n_modes), float)	
		self.EA_cluster_R = np.zeros((n_modes), float)
		self.IP_alone_R = np.zeros((n_modes), float)
		self.EA_alone_R = np.zeros((n_modes), float)
		self.P_plus_R = np.zeros((n_modes), float)
		self.P_minus_R = np.zeros((n_modes), float)	
		
		self.IP_cluster_a0 = np.zeros((n_modes), float)	
		self.EA_cluster_a0 = np.zeros((n_modes), float)
		self.IP_alone_a0 = np.zeros((n_modes), float)
		self.EA_alone_a0 = np.zeros((n_modes), float)
		self.P_plus_a0 = np.zeros((n_modes), float)
		self.P_minus_a0 = np.zeros((n_modes), float)	
		
		self.IP_cluster_a1 = np.zeros((n_modes), float)	
		self.EA_cluster_a1 = np.zeros((n_modes), float)
		self.IP_alone_a1 = np.zeros((n_modes), float)
		self.EA_alone_a1 = np.zeros((n_modes), float)
		self.P_plus_a1 = np.zeros((n_modes), float)
		self.P_minus_a1 = np.zeros((n_modes), float)	

	def __str__(self):
		tmp=''
#		for i in xrange(self.n_modes):
#			tmp += 'Mode %d\n' % (i+1)
#			for ii in xrange(self.n_displ):
#				tmp += '%12.5f %10.2f %10.2f %10.2f %12.5f %12.5f %10.2f %10.2f %10.2f %12.5f %12.5f %12.5f %12.5f\n' % (self.displ[i,ii], self.E_cluster_0[i,ii], self.E_cluster_1[i,ii], self.E_cluster_m1[i,ii], self.IP_cluster[i,ii], self.EA_cluster[i,ii], self.E_alone_0[i,ii], self.E_alone_1[i,ii], self.E_alone_m1[i,ii], self.IP_alone[i,ii], self.EA_alone[i,ii], self.P_plus[i,ii], self.P_minus[i,ii]) 

		tmp += '======================\n'
		tmp += 'Averages and variances\n'
		tmp += '======================\n'
	
		tmp += 'IP_cluster\n'
		tmp += ' Mode Freq Energy Average Std_dev DeltaE R^2 a0 a1\n'
		for i in xrange(self.n_modes):
			tmp += '%4d   %f   %e   %12.5f   %e   %e   %f   %e   %e\n' % (i+1, self.freq[i], self.energy[i], self.IP_cluster_av[i], np.sqrt(self.IP_cluster_var[i]), self.IP_cluster_dE[i], np.abs(self.IP_cluster_R[i]), self.IP_cluster_a0[i], self.IP_cluster_a1[i])
			
		tmp += 'EA_cluster\n'
		tmp += ' Mode Freq Energy Average Std_dev DeltaE R^2 a0 a1\n'
		for i in xrange(self.n_modes):
			tmp += '%4d   %f   %e   %12.5f   %e   %e   %f   %e   %e\n' % (i+1, self.freq[i], self.energy[i], self.EA_cluster_av[i], np.sqrt(self.EA_cluster_var[i]), self.EA_cluster_dE[i], np.abs(self.EA_cluster_R[i]), self.EA_cluster_a0[i], self.EA_cluster_a1[i])
			
		tmp += 'IP_alone\n'
		tmp += ' Mode Freq Energy Average Std_dev DeltaE R^2 a0 a1\n'
		for i in xrange(self.n_modes):
			tmp += '%4d   %f   %e   %12.5f   %e   %e   %f   %e   %e\n' % (i+1, self.freq[i], self.energy[i], self.IP_alone_av[i], np.sqrt(self.IP_alone_var[i]), self.IP_alone_dE[i], np.abs(self.IP_alone_R[i]), self.IP_alone_a0[i], self.IP_alone_a1[i])
			
		tmp += 'EA_alone\n'
		tmp += ' Mode Freq Energy Average Std_dev DeltaE R^2 a0 a1\n'
		for i in xrange(self.n_modes):
			tmp += '%4d   %f   %e   %12.5f   %e   %e   %f   %e   %e\n' % (i+1, self.freq[i], self.energy[i], self.EA_alone_av[i], np.sqrt(self.EA_alone_var[i]), self.EA_alone_dE[i], np.abs(self.EA_alone_R[i]), self.EA_alone_a0[i], self.EA_alone_a1[i])

		tmp += 'P_plus\n'
		tmp += ' Mode Freq Energy Average Std_dev DeltaE R^2 a0 a1\n'
		for i in xrange(self.n_modes):
			tmp += '%4d   %f   %e   %12.5f   %e   %e   %f   %e   %e\n' % (i+1, self.freq[i], self.energy[i], self.P_plus_av[i], np.sqrt(self.P_plus_var[i]), self.P_plus_dE[i], np.abs(self.P_plus_R[i]), self.P_plus_a0[i], self.P_plus_a1[i])
			
		tmp += 'P_minus\n'
		tmp += ' Mode Freq Energy Average Std_dev DeltaE R^2 a0 a1\n'
		for i in xrange(self.n_modes):
			tmp += '%4d   %f   %e   %12.5f   %e   %e   %f   %e   %e\n' % (i+1, self.freq[i], self.energy[i], self.P_minus_av[i], np.sqrt(self.P_minus_var[i]), self.P_minus_dE[i], np.abs(self.P_minus_R[i]), self.P_minus_a0[i], self.P_minus_a1[i])
			
		return tmp
		
	def Calculate_Av_Var_dE(self):
		self.IP_cluster_av = np.average(self.IP_cluster, axis=1)
		self.EA_cluster_av = np.average(self.EA_cluster, axis=1)
		self.IP_alone_av = np.average(self.IP_alone, axis=1)
		self.EA_alone_av = np.average(self.EA_alone, axis=1)
		self.P_plus_av = np.average(self.P_plus, axis=1)
		self.P_minus_av = np.average(self.P_minus, axis=1)
		
		self.IP_cluster_var = np.var(self.IP_cluster, axis=1)
		self.EA_cluster_var = np.var(self.EA_cluster, axis=1)
		self.IP_alone_var = np.var(self.IP_alone, axis=1)
		self.EA_alone_var = np.var(self.EA_alone, axis=1)
		self.P_plus_var = np.var(self.P_plus, axis=1)
		self.P_minus_var = np.var(self.P_minus, axis=1)
		
		self.IP_cluster_dE = np.max(self.IP_cluster, axis=1) - np.min(self.IP_cluster, axis=1)
		self.EA_cluster_dE = np.max(self.EA_cluster, axis=1) - np.min(self.EA_cluster, axis=1) 
		self.IP_alone_dE = np.max(self.IP_alone, axis=1) - np.min(self.IP_alone, axis=1)
		self.EA_alone_dE = np.max(self.EA_alone, axis=1) - np.min(self.EA_alone, axis=1) 
		self.P_plus_dE = np.max(self.P_plus, axis=1) - np.min(self.P_plus, axis=1) 
		self.P_minus_dE = np.max(self.P_minus, axis=1) - np.min(self.P_minus, axis=1) 
	
	def Calculate_R(self):
		S_X = np.sqrt(np.average(np.power(self.displ, 2), axis=1) - np.power(np.average(self.displ, axis=1), 2))
		
		S_Y = np.sqrt(np.average(np.power(self.IP_cluster, 2), axis=1) - np.power(self.IP_cluster_av, 2))
		S_XY = np.average(self.displ*self.IP_cluster, axis=1) - np.average(self.displ, axis=1)*np.average(self.IP_cluster, axis=1)
		self.IP_cluster_R = S_XY/(S_X*S_Y)
		
		S_Y = np.sqrt(np.average(np.power(self.EA_cluster, 2), axis=1) - np.power(self.EA_cluster_av, 2))
		S_XY = np.average(self.displ*self.EA_cluster, axis=1) - np.average(self.displ, axis=1)*np.average(self.EA_cluster, axis=1)
		self.EA_cluster_R = S_XY/(S_X*S_Y)
		
		S_Y = np.sqrt(np.average(np.power(self.IP_alone, 2), axis=1) - np.power(self.IP_alone_av, 2))
		S_XY = np.average(self.displ*self.IP_alone, axis=1) - np.average(self.displ, axis=1)*np.average(self.IP_alone, axis=1)
		self.IP_alone_R = S_XY/(S_X*S_Y)
		
		S_Y = np.sqrt(np.average(np.power(self.EA_alone, 2), axis=1) - np.power(self.EA_alone_av, 2))
		S_XY = np.average(self.displ*self.EA_alone, axis=1) - np.average(self.displ, axis=1)*np.average(self.EA_alone, axis=1)
		self.EA_alone_R = S_XY/(S_X*S_Y)
		
		S_Y = np.sqrt(np.average(np.power(self.P_plus, 2), axis=1) - np.power(self.P_plus_av, 2))
		S_XY = np.average(self.displ*self.P_plus, axis=1) - np.average(self.displ, axis=1)*np.average(self.P_plus, axis=1)
		self.P_plus_R = S_XY/(S_X*S_Y)
		
		S_Y = np.sqrt(np.average(np.power(self.P_minus, 2), axis=1) - np.power(self.P_minus_av, 2))
		S_XY = np.average(self.displ*self.P_minus, axis=1) - np.average(self.displ, axis=1)*np.average(self.P_minus, axis=1)		
		self.P_minus_R = np.power(S_XY/(S_X*S_Y),2)

class Sigma_evol(object):
	""" Class for the VBHF data"""
	def __init__(self, n_temp=1):
		self.n_temp = int(n_temp)
		
		self.IP_cluster_qu = np.zeros((n_temp), float)	
		self.EA_cluster_qu = np.zeros((n_temp), float)
		self.IP_alone_qu = np.zeros((n_temp), float)
		self.EA_alone_qu = np.zeros((n_temp), float)
		self.P_plus_qu = np.zeros((n_temp), float)
		self.P_minus_qu = np.zeros((n_temp), float)
		
		self.IP_cluster_cl = np.zeros((n_temp), float)	
		self.EA_cluster_cl = np.zeros((n_temp), float)
		self.IP_alone_cl = np.zeros((n_temp), float)
		self.EA_alone_cl = np.zeros((n_temp), float)
		self.P_plus_cl = np.zeros((n_temp), float)
		self.P_minus_cl = np.zeros((n_temp), float)
		
	def __str__(self):
		tmp = ''
		
		tmp += "IP_cluster\n"
		for i in xrange(self.n_temp):
			tmp += "%d %e %e\n" % (i+1, self.IP_cluster_qu[i], self.IP_cluster_cl[i])

		tmp += "EA_cluster\n"
		for i in xrange(self.n_temp):
			tmp += "%d %e %e\n" % (i+1, self.EA_cluster_qu[i], self.EA_cluster_cl[i])
			
		tmp += "IP_alone\n"
		for i in xrange(self.n_temp):
			tmp += "%d %e %e\n" % (i+1, self.IP_alone_qu[i], self.IP_alone_cl[i])
			
		tmp += "EA_alone\n"
		for i in xrange(self.n_temp):
			tmp += "%d %e %e\n" % (i+1, self.EA_alone_qu[i], self.EA_alone_cl[i])	
				
		tmp += "P_plus\n"
		for i in xrange(self.n_temp):
			tmp += "%d %e %e\n" % (i+1, self.P_plus_qu[i], self.P_plus_cl[i])
			
		tmp += "P_minus\n"
		for i in xrange(self.n_temp):
			tmp += "%d %e %e\n" % (i+1, self.P_minus_qu[i], self.P_minus_cl[i])
			
		return tmp	
			
	def Calculate_Sigma(self, data):
		for temperature in xrange(self.n_temp):
			tmp = 1/(sp.tanh(data.energy/(2*K_B*(temperature+1))))
			self.IP_cluster_qu[temperature] = np.sqrt(np.sum((data.IP_cluster_a1*data.IP_cluster_a1)/2 * tmp))
			self.EA_cluster_qu[temperature] = np.sqrt(np.sum((data.EA_cluster_a1*data.EA_cluster_a1)/2 * tmp))
			self.IP_alone_qu[temperature] = np.sqrt(np.sum((data.IP_alone_a1*data.IP_alone_a1)/2 * tmp))
			self.EA_alone_qu[temperature] = np.sqrt(np.sum((data.EA_alone_a1*data.EA_alone_a1)/2 * tmp))
			self.P_plus_qu[temperature] = np.sqrt(np.sum((data.P_plus_a1*data.P_plus_a1)/2 * tmp))
			self.P_minus_qu[temperature] = np.sqrt(np.sum((data.P_minus_a1*data.P_minus_a1)/2 * tmp))
			
			tmp = 2*K_B*(temperature+1)/data.energy
			self.IP_cluster_cl[temperature] = np.sqrt(np.sum((data.IP_cluster_a1*data.IP_cluster_a1)/2 * tmp))
			self.EA_cluster_cl[temperature] = np.sqrt(np.sum((data.EA_cluster_a1*data.EA_cluster_a1)/2 * tmp))
			self.IP_alone_cl[temperature] = np.sqrt(np.sum((data.IP_alone_a1*data.IP_alone_a1)/2 * tmp))
			self.EA_alone_cl[temperature] = np.sqrt(np.sum((data.EA_alone_a1*data.EA_alone_a1)/2 * tmp))
			self.P_plus_cl[temperature] = np.sqrt(np.sum((data.P_plus_a1*data.P_plus_a1)/2 * tmp))
			self.P_minus_cl[temperature] = np.sqrt(np.sum((data.P_minus_a1*data.P_minus_a1)/2 * tmp))
			
			"""tmp = 0.0
			tmp2 = 0.0
			for i in xrange(144):
				tmp += (2*K_B*(temperature+1)/data.energy[i])*(data.P_plus_a1[i]*data.P_plus_a1[i])/2
				tmp2 += (data.P_plus_a1[i]*data.P_plus_a1[i])/2
			print temperature, np.sqrt(tmp)
			
			self.P_plus_cl[temperature] = np.sqrt(tmp)"""
			
			
			
#=====================================================================
#--------------------------Data Read/Write----------------------------
#=====================================================================

def Read_Data(data):
	""" Read VBHF data file """
	input_file = "%s.dat" % (FILE)
	
	try:
		finput = open(input_file, 'r')
	except:
		print "Could not open %s" % (input_file)
		sys.exit(1)

	for i in xrange(data.n_modes):
			line = finput.readline()
			words = line.split()
			while (len(words)==0):
				line = finput.readline()
				words = line.split()
						
			if (len(words)!=0 and (words[0]=="Mode" or words[0]=="mode")):
				for ii in xrange(data.n_displ):
					line = finput.readline()
					words = line.split()
					while (len(words)==0):
						line = finput.readline()
						words = line.split()
					
					data.displ[i,ii] = float(words[0])
					data.E_cluster_0[i,ii] = float(words[1])
					data.E_cluster_1[i,ii] = float(words[2])
					data.E_cluster_m1[i,ii] = float(words[3])

					data.E_alone_0[i,ii] = float(words[6])
					data.E_alone_1[i,ii] = float(words[7])
					data.E_alone_m1[i,ii] = float(words[8])

					data.IP_cluster[i,ii] = float(words[4])
					data.EA_cluster[i,ii] = float(words[5])
					data.IP_alone[i,ii] = float(words[9])
					data.EA_alone[i,ii] = float(words[10])
					data.P_plus[i,ii] = float(words[11])
					data.P_minus[i,ii] = float(words[12])
				
	finput.close()

	input_file = "%s.freq" % (FILE)

	try:
		finput = open(input_file, 'r')
	except:
		print "Could not open %s" % (input_file)
		sys.exit(1)
		
	for i in xrange(data.n_modes):
			line = finput.readline()
			words = line.split()
			while (len(words)==0):
				line = finput.readline()
				words = line.split()
						
			if len(words)==3 :
				data.freq[i] = float(words[1])
				data.energy[i] = float(words[2])
				
	finput.close()
#	print "Reading of input file done!"

def Read_Fit(data):
	""" Read fit data files """
	
	list=["data_fit.IP_cluster", "data_fit.EA_cluster", \
	"data_fit.IP_alone", "data_fit.EA_alone", \
	"data_fit.P_minus", "data_fit.P_plus"]
	
	for input_file in list:
		input = True
		try:
			finput = open(input_file, 'r')
		except:
			input = False
		
		if input:
			for i in xrange(data.n_modes):
					line = finput.readline()
					words = line.split()
					while (len(words)==0):
						line = finput.readline()
						words = line.split()
								
					if (len(words)==3):
						if input_file=="data_fit.P_minus" and abs(data.P_minus_R[i]) >= R_MIN:
							data.P_minus_a0[i] = float(words[1])
							data.P_minus_a1[i] = float(words[2])
							
						elif input_file=="data_fit.P_plus" and abs(data.P_plus_R[i]) >= R_MIN:
							data.P_plus_a0[i] = float(words[1])
							data.P_plus_a1[i] = float(words[2])

						elif input_file=="data_fit.IP_cluster" and abs(data.IP_cluster_R[i]) >= R_MIN:
							data.IP_cluster_a0[i] = float(words[1])
							data.IP_cluster_a1[i] = float(words[2])
							
						elif input_file=="data_fit.EA_cluster" and abs(data.EA_cluster_R[i]) >= R_MIN:
							data.EA_cluster_a0[i] = float(words[1])
							data.EA_cluster_a1[i] = float(words[2])
							
						elif input_file=="data_fit.IP_alone" and abs(data.IP_alone_R[i]) >= R_MIN:
							data.IP_alone_a0[i] = float(words[1])
							data.IP_alone_a1[i] = float(words[2])
							
						elif input_file=="data_fit.EA_alone" and abs(data.EA_alone_R[i]) >= R_MIN:
							data.EA_alone_a0[i] = float(words[1])
							data.EA_alone_a1[i] = float(words[2])

			finput.close()

def Write_Gnuplot_dat(data):
	try :
		os.mkdir("plots")
	except:
		pass
	
	for i in xrange(data.n_modes):
		output_file = "plots/mode_%d.dat" % (i+1)
		try:
			foutput = open(output_file, 'w')
		except:
			print "Could not open %s" % (output_file)
			exit(1)
		
		tmp = ''
		for ii in xrange(data.n_displ):
			tmp += '%12.5e %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n' % (data.displ[i,ii], data.IP_cluster[i,ii], data.EA_cluster[i,ii], data.IP_alone[i,ii], data.EA_alone[i,ii], data.P_plus[i,ii], data.P_minus[i,ii]) 
		
		foutput.write(tmp)
		foutput.close()

def Write_Gnuplot_plt(data, type):
	try :
		os.mkdir("plots/IP_cluster")
		os.mkdir("plots/EA_cluster")
		os.mkdir("plots/IP_alone")
		os.mkdir("plots/EA_alone")
		os.mkdir("plots/P_plus")
		os.mkdir("plots/P_minus")
	except:
		pass
	
	for i in xrange(data.n_modes):
		output_file = "plots/mode_%d.plt" % (i+1)
		try:
			foutput = open(output_file, 'w')
		except:
			print "Could not open %s" % (output_file)
			sys.exit(1)
		
		if (type=="epslatex"):
			tmp = 'reset\n'
			tmp += 'set xlabel \'Displacement (normal coordinates)\'\n'
			tmp += 'set ylabel \'$E\\ (eV)$\'\n'
			tmp += 'set terminal epslatex size 5,3\n\n'
			tmp += '##################################\n\n'
			
			tmp += 'set output \'IP_cluster/mode_%d.eps\'\n' % (i+1)
			tmp += 'f(x) = a*x+b\n'
			tmp += 'fit f(x) \'mode_%d.dat\' using 1:2 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:2 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'

			tmp += 'set output \'EA_cluster/mode_%d.eps\'\n' % (i+1)
			tmp += 'fit f(x) \'mode_%d.dat\' using 1:3 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:3 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'
			
			tmp += 'set output \'IP_alone/mode_%d.eps\'\n' % (i+1)
			tmp += 'fit f(x) \'mode_%d.dat\' using 1:4 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:4 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'
			
			tmp += 'set output \'EA_alone/mode_%d.eps\'\n' % (i+1)
			tmp += 'fit f(x) \'mode_%d.dat\' using 1:5 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:5 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'

			tmp += 'set output \'P_plus/mode_%d.eps\'\n' % (i+1)
			tmp += 'fit f(x) \'mode_%d.dat\' using 1:6 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:6 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'
			
			tmp += 'set output \'P_minus/mode_%d.eps\'\n' % (i+1)
			tmp += 'fit f(x) \'mode_%d.dat\' using 1:7 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:7 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'

		if (type=="png"):
			tmp = 'reset\n'
			tmp += 'set xlabel \'Displacement (normal coordinates)\'\n'
			tmp += 'set ylabel \'E (eV)\'\n'
			tmp += 'set terminal png size 640,480 font \'/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif.ttf\'\n\n'
			tmp += '##################################\n\n'
			
			tmp += 'set output \'IP_cluster/mode_%d.png\'\n' % (i+1)
			tmp += 'f(x) = a*x+b\n'
			tmp += 'fit f(x) \'mode_%d.dat\' using 1:2 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:2 lw 6 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'

			tmp += 'set output \'EA_cluster/mode_%d.png\'\n' % (i+1)
			tmp += 'fit f(x) \'mode_%d.dat\' using 1:3 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:3 lw 6 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'
			
			tmp += 'set output \'IP_alone/mode_%d.png\'\n' % (i+1)
			tmp += 'fit f(x) \'mode_%d.dat\' using 1:4 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:4 lw 6 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'
			
			tmp += 'set output \'EA_alone/mode_%d.png\'\n' % (i+1)
			tmp += 'fit f(x) \'mode_%d.dat\' using 1:5 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:5 lw 6 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'

			tmp += 'set output \'P_plus/mode_%d.png\'\n' % (i+1)
			tmp += 'fit f(x) \'mode_%d.dat\' using 1:6 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:6 lw 6 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'
			
			tmp += 'set output \'P_minus/mode_%d.png\'\n' % (i+1)
			tmp += 'fit f(x) \'mode_%d.dat\' using 1:7 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:7 lw 6 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'
		
		foutput.write(tmp)
		foutput.close()

if __name__ == '__main__':
	""" 
		This analyze a file containing the polarization energies, IP, EA, etc.
		Then, it: 
			- creates Gnuplot files for studying the fluctuations of the energies 
		  	  along each mode 
			- calcultes averages, variances and energy differences for each normal
		  	  mode
			- calculates the linearity coefficient R
			- run the Gnuplot scripts, and collect the fit values
			
		Finally, it calculates the evolution of sigma for given temperatures at
		the classical and quantum level.
		
	"""
	t1 = time.clock()
	
	(n_modes, n_displ) = (144, 11)
	
	# Read the file containing polarization energies and the file
	# containing the frequencies and energies of the modes
	qs=VBHF_Data(n_modes, n_displ)
	Read_Data(qs)
	
	# Calculate averages, variances and deltaE for information
	qs.Calculate_Av_Var_dE()
	
	# Calculate R coefficients
	qs.Calculate_R()
	
	# Create Gnuplot input files
	Write_Gnuplot_dat(qs)
	Write_Gnuplot_plt(qs, "png")

	# Copy and run the Gnuplot scripts to calculate the fits	
	try:
		shutil.copy("src%sQS_collect_fit_parameters_red.sh" % (os.sep), "plots")
		os.chdir("plots")
		os.system("chmod +x QS_collect_fit_parameters_red.sh")
	except:
		print "[ERROR] Problem when copying file QS_collect_fit_parameters_red.sh"
		sys.exit(1)
	os.system("./QS_collect_fit_parameters_red.sh")
	
	# Read the results of the fits
	Read_Fit(qs)

	# Calculate the evolution of sigma
	qs_sig = Sigma_evol(N_TEMP)
	qs_sig.Calculate_Sigma(qs)
	
	"""output_file = "%s.an" % (FILE)
	try:
		foutput = open(output_file, 'w')
	except:
		print "Could not open %s" % (output_file)
		sys.exit(1)
	tmp = print qs
	foutput.write(tmp)
	foutput.close()
	
	output_file = "%s.sigma" % (FILE)
	try:
		foutput = open(output_file, 'w')
	except:
		print "Could not open %s" % (output_file)
		sys.exit(1)
	tmp = print qs_sig
	foutput.write(tmp)
	foutput.close()"""
	
	print qs
	print qs_sig
	
	t2 = time.clock()
	print t2-t1
