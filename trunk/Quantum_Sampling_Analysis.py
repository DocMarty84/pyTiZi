#!/usr/bin/env python

import string
import os, shutil, sys
import math
import numpy
import time
import copy

R_MIN = 0.97

#=====================================================================
#-------------------------Class Definition----------------------------
#=====================================================================

class VBHF_Data(object):
	""" Class for the VBHF data"""
	def __init__(self, n_modes=1, n_displ=1):
		self.n_modes = n_modes	# Number of vibrational modes
		self.n_displ = n_displ	# Number of displacements, including equilibrium geometry
		
		self.displ = numpy.zeros((n_modes, n_displ), float)			# Displacement in normal modes coordinates
		
		self.E_cluster_0 = numpy.zeros((n_modes, n_displ), float)	# VBHF energy for neutral cluster
		self.E_cluster_1 = numpy.zeros((n_modes, n_displ), float)	# VBHF energy for cluster with charge=1
		self.E_cluster_m1 = numpy.zeros((n_modes, n_displ), float)	# VBHF energy for cluster with charge=-1

		self.E_alone_0 = numpy.zeros((n_modes, n_displ), float)		# VBHF energy for alone neutral molecule
		self.E_alone_1 = numpy.zeros((n_modes, n_displ), float)		# VBHF energy for alone molecule with charge=1
		self.E_alone_m1 = numpy.zeros((n_modes, n_displ), float)	# VBHF energy for alone molecule with charge=-1

		self.IP_cluster = numpy.zeros((n_modes, n_displ), float)	# VBHF energy for alone molecule with charge=1
		self.EA_cluster = numpy.zeros((n_modes, n_displ), float)	# VBHF energy for alone molecule with charge=-1
		self.IP_alone = numpy.zeros((n_modes, n_displ), float)		# VBHF energy for alone molecule with charge=1
		self.EA_alone = numpy.zeros((n_modes, n_displ), float)		# VBHF energy for alone molecule with charge=-1
		self.P_plus = numpy.zeros((n_modes, n_displ), float)		# VBHF energy for alone molecule with charge=1
		self.P_minus = numpy.zeros((n_modes, n_displ), float)		# VBHF energy for alone molecule with charge=-1	

		self.IP_cluster_av = numpy.zeros((n_modes), float)	
		self.EA_cluster_av = numpy.zeros((n_modes), float)
		self.IP_alone_av = numpy.zeros((n_modes), float)
		self.EA_alone_av = numpy.zeros((n_modes), float)
		self.P_plus_av = numpy.zeros((n_modes), float)
		self.P_minus_av = numpy.zeros((n_modes), float)	
		
		self.IP_cluster_var = numpy.zeros((n_modes), float)	
		self.EA_cluster_var = numpy.zeros((n_modes), float)
		self.IP_alone_var = numpy.zeros((n_modes), float)
		self.EA_alone_var = numpy.zeros((n_modes), float)
		self.P_plus_var = numpy.zeros((n_modes), float)
		self.P_minus_var = numpy.zeros((n_modes), float)
	
		self.IP_cluster_dE = numpy.zeros((n_modes), float)	
		self.EA_cluster_dE = numpy.zeros((n_modes), float)
		self.IP_alone_dE = numpy.zeros((n_modes), float)
		self.EA_alone_dE = numpy.zeros((n_modes), float)
		self.P_plus_dE = numpy.zeros((n_modes), float)
		self.P_minus_dE = numpy.zeros((n_modes), float)	
		
		self.IP_cluster_R = numpy.zeros((n_modes), float)	
		self.EA_cluster_R = numpy.zeros((n_modes), float)
		self.IP_alone_R = numpy.zeros((n_modes), float)
		self.EA_alone_R = numpy.zeros((n_modes), float)
		self.P_plus_R = numpy.zeros((n_modes), float)
		self.P_minus_R = numpy.zeros((n_modes), float)	
		
		self.IP_cluster_a0 = numpy.zeros((n_modes), float)	
		self.EA_cluster_a0 = numpy.zeros((n_modes), float)
		self.IP_alone_a0 = numpy.zeros((n_modes), float)
		self.EA_alone_a0 = numpy.zeros((n_modes), float)
		self.P_plus_a0 = numpy.zeros((n_modes), float)
		self.P_minus_a0 = numpy.zeros((n_modes), float)	
		
		self.IP_cluster_a1 = numpy.zeros((n_modes), float)	
		self.EA_cluster_a1 = numpy.zeros((n_modes), float)
		self.IP_alone_a1 = numpy.zeros((n_modes), float)
		self.EA_alone_a1 = numpy.zeros((n_modes), float)
		self.P_plus_a1 = numpy.zeros((n_modes), float)
		self.P_minus_a1 = numpy.zeros((n_modes), float)	

	def __str__(self):
		tmp=''
#		for i in xrange(self.n_modes):
#			tmp += 'Mode %d\n' % (i+1)
#			for ii in xrange(self.n_displ):
#				tmp += '%6.2f %10.2f %10.2f %10.2f %12.5f %12.5f %10.2f %10.2f %10.2f %12.5f %12.5f %12.5f %12.5f\n' % (self.displ[i,ii], self.E_cluster_0[i,ii], self.E_cluster_1[i,ii], self.E_cluster_m1[i,ii], self.IP_cluster[i,ii], self.EA_cluster[i,ii], self.E_alone_0[i,ii], self.E_alone_1[i,ii], self.E_alone_m1[i,ii], self.IP_alone[i,ii], self.EA_alone[i,ii], self.P_plus[i,ii], self.P_minus[i,ii]) 

		tmp += '======================\n'
		tmp += 'Averages and variances\n'
		tmp += '======================\n'
	
		tmp += 'IP_cluster\n'
		tmp += 'Mode Average Std_dev DeltaE R a0 a1\n'
		for i in xrange(self.n_modes):
			tmp += '%4d %12.5f   %e   %e   %f   %e   %e\n' % (i+1, self.IP_cluster_av[i], numpy.sqrt(self.IP_cluster_var[i]), self.IP_cluster_dE[i], numpy.abs(self.IP_cluster_R[i]), self.IP_cluster_a0[i], self.IP_cluster_a1[i])
			
		tmp += 'EA_cluster\n'
		tmp += 'Mode Average Std_dev DeltaE R a0 a1\n'
		for i in xrange(self.n_modes):
			tmp += '%4d %12.5f   %e   %e   %f   %e   %e\n' % (i+1, self.EA_cluster_av[i], numpy.sqrt(self.EA_cluster_var[i]), self.EA_cluster_dE[i], numpy.abs(self.EA_cluster_R[i]), self.EA_cluster_a0[i], self.EA_cluster_a1[i])
			
		tmp += 'IP_alone\n'
		tmp += 'Mode Average Std_dev DeltaE R a0 a1\n'
		for i in xrange(self.n_modes):
			tmp += '%4d %12.5f   %e   %e   %f   %e   %e\n' % (i+1, self.IP_alone_av[i], numpy.sqrt(self.IP_alone_var[i]), self.IP_alone_dE[i], numpy.abs(self.IP_alone_R[i]), self.IP_alone_a0[i], self.IP_alone_a0[i])
			
		tmp += 'EA_alone\n'
		tmp += 'Mode Average Std_dev DeltaE R a0 a1\n'
		for i in xrange(self.n_modes):
			tmp += '%4d %12.5f   %e   %e   %f   %e   %e\n' % (i+1, self.EA_alone_av[i], numpy.sqrt(self.EA_alone_var[i]), self.EA_alone_dE[i], numpy.abs(self.EA_alone_R[i]), self.EA_alone_a0[i], self.EA_alone_a0[i])

		tmp += 'P_plus\n'
		tmp += 'Mode Average Std_dev DeltaE R a0 a1\n'
		for i in xrange(self.n_modes):
			tmp += '%4d %12.5f   %e   %e   %f   %e   %e\n' % (i+1, self.P_plus_av[i], numpy.sqrt(self.P_plus_var[i]), self.P_plus_dE[i], numpy.abs(self.P_plus_R[i]), self.P_plus_a0[i], self.P_plus_a0[i])
			
		tmp += 'P_minus\n'
		tmp += 'Mode Average Std_dev DeltaE R a0 a1\n'
		for i in xrange(self.n_modes):
			tmp += '%4d %12.5f   %e   %e   %f   %e   %e\n' % (i+1, self.P_minus_av[i], numpy.sqrt(self.P_minus_var[i]), self.P_minus_dE[i], numpy.abs(self.P_minus_R[i]), self.P_minus_a0[i], self.P_minus_a0[i])
			
		return tmp
		
	def Calculate_Av_Var_dE(self):
		self.IP_cluster_av = numpy.average(self.IP_cluster, axis=1)
		self.EA_cluster_av = numpy.average(self.EA_cluster, axis=1)
		self.IP_alone_av = numpy.average(self.IP_alone, axis=1)
		self.EA_alone_av = numpy.average(self.EA_alone, axis=1)
		self.P_plus_av = numpy.average(self.P_plus, axis=1)
		self.P_minus_av = numpy.average(self.P_minus, axis=1)
		
		self.IP_cluster_var = numpy.var(self.IP_cluster, axis=1)
		self.EA_cluster_var = numpy.var(self.EA_cluster, axis=1)
		self.IP_alone_var = numpy.var(self.IP_alone, axis=1)
		self.EA_alone_var = numpy.var(self.EA_alone, axis=1)
		self.P_plus_var = numpy.var(self.P_plus, axis=1)
		self.P_minus_var = numpy.var(self.P_minus, axis=1)
		
		self.IP_cluster_dE = numpy.max(self.IP_cluster, axis=1) - numpy.min(self.IP_cluster, axis=1)
		self.EA_cluster_dE = numpy.max(self.EA_cluster, axis=1) - numpy.min(self.EA_cluster, axis=1) 
		self.IP_alone_dE = numpy.max(self.IP_alone, axis=1) - numpy.min(self.IP_alone, axis=1)
		self.EA_alone_dE = numpy.max(self.EA_alone, axis=1) - numpy.min(self.EA_alone, axis=1) 
		self.P_plus_dE = numpy.max(self.P_plus, axis=1) - numpy.min(self.P_plus, axis=1) 
		self.P_minus_dE = numpy.max(self.P_minus, axis=1) - numpy.min(self.P_minus, axis=1) 
	
	def Calculate_R(self):
		S_X = numpy.sqrt(numpy.average(numpy.power(self.displ, 2), axis=1) - numpy.power(numpy.average(self.displ, axis=1), 2))
		
		S_Y = numpy.sqrt(numpy.average(numpy.power(self.IP_cluster, 2), axis=1) - numpy.power(self.IP_cluster_av, 2))
		S_XY = numpy.average(self.displ*self.IP_cluster, axis=1) - numpy.average(self.displ, axis=1)*numpy.average(self.IP_cluster, axis=1)
		self.IP_cluster_R = S_XY/(S_X*S_Y)
		
		S_Y = numpy.sqrt(numpy.average(numpy.power(self.EA_cluster, 2), axis=1) - numpy.power(self.EA_cluster_av, 2))
		S_XY = numpy.average(self.displ*self.EA_cluster, axis=1) - numpy.average(self.displ, axis=1)*numpy.average(self.EA_cluster, axis=1)
		self.EA_cluster_R = S_XY/(S_X*S_Y)
		
		S_Y = numpy.sqrt(numpy.average(numpy.power(self.IP_alone, 2), axis=1) - numpy.power(self.IP_alone_av, 2))
		S_XY = numpy.average(self.displ*self.IP_alone, axis=1) - numpy.average(self.displ, axis=1)*numpy.average(self.IP_alone, axis=1)
		self.IP_alone_R = S_XY/(S_X*S_Y)
		
		S_Y = numpy.sqrt(numpy.average(numpy.power(self.EA_alone, 2), axis=1) - numpy.power(self.EA_alone_av, 2))
		S_XY = numpy.average(self.displ*self.EA_alone, axis=1) - numpy.average(self.displ, axis=1)*numpy.average(self.EA_alone, axis=1)
		self.EA_alone_R = S_XY/(S_X*S_Y)
		
		S_Y = numpy.sqrt(numpy.average(numpy.power(self.P_plus, 2), axis=1) - numpy.power(self.P_plus_av, 2))
		S_XY = numpy.average(self.displ*self.P_plus, axis=1) - numpy.average(self.displ, axis=1)*numpy.average(self.P_plus, axis=1)
		self.P_plus_R = S_XY/(S_X*S_Y)
		
		S_Y = numpy.sqrt(numpy.average(numpy.power(self.P_minus, 2), axis=1) - numpy.power(self.P_minus_av, 2))
		S_XY = numpy.average(self.displ*self.P_minus, axis=1) - numpy.average(self.displ, axis=1)*numpy.average(self.P_minus, axis=1)		
		self.P_minus_R = S_XY/(S_X*S_Y)
		
#=====================================================================
#--------------------------Data Read/Write----------------------------
#=====================================================================

def Read_Data(input_file, data):
	""" Read VBHF data file """
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
#	print "Reading of input file done!"

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
			exit(1)
		
		if (type=="epslatex"):
			tmp = 'reset\n'
			tmp += 'set xlabel \'Displacement (normal coordinates)\'\n'
			tmp += 'set ylabel \'$E\\ (eV)$\'\n'
			tmp += 'set terminal epslatex size 5,3\n\n'
			tmp += '##################################\n\n'
			
			tmp += 'set output \'IP_cluster/mode_%d.eps\'\n' % (i+1)
			tmp += 'f(x) = a*x+b\n'
			tmp += '#fit f(x) \'mode_%d.dat\' using 1:2 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:2 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'

			tmp += 'set output \'EA_cluster/mode_%d.eps\'\n' % (i+1)
			tmp += '#fit f(x) \'mode_%d.dat\' using 1:3 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:3 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'
			
			tmp += 'set output \'IP_alone/mode_%d.eps\'\n' % (i+1)
			tmp += '#fit f(x) \'mode_%d.dat\' using 1:4 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:4 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'
			
			tmp += 'set output \'EA_alone/mode_%d.eps\'\n' % (i+1)
			tmp += '#fit f(x) \'mode_%d.dat\' using 1:5 via a,b\n' % (i+1)
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
			tmp += '#fit f(x) \'mode_%d.dat\' using 1:2 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:2 lw 6 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'

			tmp += 'set output \'EA_cluster/mode_%d.png\'\n' % (i+1)
			tmp += '#fit f(x) \'mode_%d.dat\' using 1:3 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:3 lw 6 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'
			
			tmp += 'set output \'IP_alone/mode_%d.png\'\n' % (i+1)
			tmp += '#fit f(x) \'mode_%d.dat\' using 1:4 via a,b\n' % (i+1)
			tmp += 'plot \'mode_%d.dat\' using 1:4 lw 6 title \'Mode %d\'\n' % (i+1, i+1)
			tmp += 'set output\n\n'
			
			tmp += 'set output \'EA_alone/mode_%d.png\'\n' % (i+1)
			tmp += '#fit f(x) \'mode_%d.dat\' using 1:5 via a,b\n' % (i+1)
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
			
		The coefficient of the linear fit will be obtained while running the
		Gnuplot scripts. They can be collected with the script 
		collect_fit_parameters.sh.
		
	"""
	t1 = time.clock()
	
	(n_modes, n_displ) = (144, 11)
	
	qs=VBHF_Data(n_modes, n_displ)
	Read_Data("AM1.dat", qs)
	
	Write_Gnuplot_dat(qs)
	Write_Gnuplot_plt(qs, "png")

	qs.Calculate_Av_Var_dE()
	qs.Calculate_R()
	
	try:
		shutil.copy("src%scollect_fit_parameters_red.sh" % (os.sep), "plots")
		os.chdir("plots")
		os.system("chmod +x collect_fit_parameters_red.sh" % (os.sep))
	except:
		print "[ERROR] Problem when copying file collect_fit_parameters_red.sh"
	
	os.system("./collect_fit_parameters_red.sh")
	
	Read_Fit(qs)
	
	print qs
	
	t2 = time.clock()
	print t2-t1
