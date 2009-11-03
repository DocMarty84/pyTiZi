#!/usr/bin/env python

import string
import os, shutil, sys
import time
import copy

sys.path.append("lib")
import input_reading_MD, molecular_system, list_manipulation, cluster_file_creation

def ls(path, ext):
	return [fichier for fichier in os.listdir(path)
								if os.path.splitext(fichier)[1] == ext]
		
def FileTypeCheck(file_type):
	if file_type == "tinker":
		return ("tinker", ".arc")
	else:
		print "[ERROR] File type not recognized. Aborting..."
		sys.exit(1)

def ClearScreen(numlines=100):
	"""Clear the console.
	numlines is an optional argument used only as a fall-back."""

	if os.name == "posix": # Unix/Linux/MacOS/BSD/etc
		os.system('clear')
	elif os.name in ("nt", "dos", "ce"): # DOS/Windows
		os.system('CLS')
	else: # Fallback for other operating systems
		print '\n' * numlines

			
if __name__ == '__main__':
	t1 = time.clock()

	input_instructions = "test_lyra"
	
	##############################################
	#                                            #
	#   Verbosity of the program:                #
	#       0: nothing is displayed              #
	#       1: errors                            #
	#       2: warning (default)                 #
	#       3: informations                      #
	#       4: additional informations           #
	#          (recommended for large systems)   #
	#                                            #
	##############################################
	
	verb = 3

	# ==========================
	# Import Psyco if available
	# ==========================
	try:
		import psyco
		psyco.profile()
		if verb > 3:
			print "[INFO] Psyco works!"
	except ImportError:
		if verb > 3:
			print "[INFO] Psyco not found! Visit http://psyco.sourceforge.net/"
	ClearScreen()
	
	# ==========================
	# Reading instructions file
	# ==========================
	try:
		finstructions = open(input_instructions, 'r')
	except:
		if verb > 0:
			print "[ERROR] Could not find %s file. Aborting.\n" % (input_instructions)
		sys.exit(1)
	
	project_name_ok = False
	input_dir_ok = False
	file_type_ok = False
	input_dir_local_ok = False
	molecules_to_analyze_ok = False
	pbc_ok = False
	molecules_for_J_ok = False
	cutoff_ok = False
	username_cluster_ok = False
	dir_cluster_ok = False
	input_dir_cluster_ok = False
	output_dir_cluster_ok = False
	scratch_dir_cluster_ok = False
	location_cluster_ok = False
	zindo_dir_cluster_ok = False
	
	for line in finstructions:
		words = line.split()
		if len(words) != 0:
			words[0] = words[0].strip(string.punctuation + string.whitespace)
			first_word = string.upper(words[0])
			
			if first_word=="PROJECT_NAME":
				project_name = words[1]
				project_name_ok = True
			
			elif first_word=="INPUT_DIR":
				input_dir = words[1]
				input_dir_ok = True
				
			elif first_word=="FILE_TYPE":
				file_type = string.lower(words[1])
				file_type_ok = True
			
			elif first_word=="MOLECULES_TO_ANALYZE":
	#			cut = re.compile('[;, ]')
	#			molecules_to_analyze = [ cut.split(x) for x in words[1:] ]
				a = [ x.split(",") for x in words[1:] ]
				molecules_to_analyze = [ string.lower(item) for sublist in a for item in sublist ]
				molecules_to_analyze_ok = True
			
			elif first_word=="PBC":
				a = [ x.split(",") for x in words[1:] ]
				pbc = [ string.lower(item) for sublist in a for item in sublist ]
				pbc_ok = True
			
			elif first_word=="MOLECULES_FOR_J":
				a = [ x.split(",") for x in words[1:] ]
				molecules_for_J = [ string.lower(item) for sublist in a for item in sublist ]
				molecules_for_J_ok = True
			
			elif first_word=="CUTOFF":
				cutoff = float(words[1])
				cutoff_ok = True
			
			elif first_word=="USERNAME_CLUSTER":
				username_cluster = words[1]
				username_cluster_ok = True

			elif first_word=="DIR_CLUSTER":
				dir_cluster = words[1]
				dir_cluster_ok = True

			elif first_word=="INPUT_DIR_CLUSTER":
				input_dir_cluster = words[1]
				input_dir_cluster_ok = True
				
			elif first_word=="OUTPUT_DIR_CLUSTER":
				output_dir_cluster = words[1]
				output_dir_cluster_ok = True
			
			elif first_word=="SCRATCH_DIR_CLUSTER":
				scratch_dir_cluster = words[1]
				scratch_dir_cluster_ok = True
				
			elif first_word=="LOCATION_CLUSTER":
				location_cluster = words[1]
				location_cluster_ok = True
			
			elif first_word=="ZINDO_DIR_CLUSTER":
				zindo_dir_cluster = words[1]
				zindo_dir_cluster_ok = True
				
			elif first_word=="END":
				break

	finstructions.close()

	quit = False
	
	# ==========================
	#   Variables obligatoires
	# ==========================
	
	while not project_name_ok:
		print "Variable PROJECT_NAME not specified!"
		a = raw_input("PROJECT_NAME = ")
		while len(a) == 0:
			a = raw_input("PROJECT_NAME = ")
			
		project_name = copy.copy(a)
		project_name_ok = True

	while not file_type_ok:
		print "Variable FILE_TYPE not specified!"
		a = raw_input("FILE_TYPE = ")
		while len(a) == 0:
			a = raw_input("FILE_TYPE = ")
			
		file_type = string.lower(a)
		file_type_ok = True
	
	while not molecules_to_analyze_ok:
		print "Variable MOLECULES_TO_ANALYZE not specified!"
		a = raw_input("MOLECULES_TO_ANALYZE = ")
		while len(a) == 0:
			a = raw_input("MOLECULES_TO_ANALYZE = ")
		
		a = [ x.split(",") for x in a ]
		molecules_to_analyze = [ string.lower(item) for sublist in a for item in sublist ]
		molecules_to_analyze_ok = True
		
	while not pbc_ok:
		print "Variable PBC not specified!"
		a = raw_input("PBC = ")
		while len(a) == 0:
			a = raw_input("PBC = ")

		a = [ x.split(",") for x in a ]
		pbc = [ string.lower(item) for sublist in a for item in sublist ]
		pbc_ok = True

	while not molecules_for_J_ok:
		print "Variable MOLECULES_FOR_J not specified!"
		a = raw_input("MOLECULES_FOR_J = ")
		while len(a) == 0:
			a = raw_input("MOLECULES_FOR_J = ")

		a = [ x.split(",") for x in a ]
		molecules_for_J = [ string.lower(item) for sublist in a for item in sublist ]
		molecules_for_J_ok = True

	while not cutoff_ok:
		print "Variable CUTOFF not specified!"
		a = raw_input("CUTOFF = ")
		while len(a) == 0:
			a = raw_input("CUTOFF = ")
			
		cutoff = float(a)
		cutoff_ok = True

	while not username_cluster_ok:
		print "Variable USERNAME_CLUSTER not specified!"
		a = raw_input("USERNAME_CLUSTER = ")
		while len(a) == 0:
			a = raw_input("USERNAME_CLUSTER = ")
			
		username_cluster = copy.copy(a)
		username_cluster_ok = True
		
	while not location_cluster_ok:
		print "Variable LOCATION_CLUSTER not specified!"
		a = raw_input("LOCATION_CLUSTER = ")
		while len(a) == 0:
			a = raw_input("LOCATION_CLUSTER = ")
			
		location_cluster = copy.copy(a)
		location_cluster_ok = True

	# ==========================
	# Variables non obligatoires
	# ==========================

	while not input_dir_ok:
		print "Variable INPUT_DIR not specified!"
		a = raw_input("INPUT_DIR [input_data] = ")
		if len(a) != 0:
			input_dir = copy.copy(a)
		else:
			input_dir = "input_data"
		input_dir_ok = True

	while not dir_cluster_ok:
		print "Variable DIR_CLUSTER not specified!"
		a = raw_input("DIR_CLUSTER [/output/%s/temp/%s] = " % (username_cluster, project_name))
		if len(a) != 0:
			dir_cluster = copy.copy(a)
		else:
			dir_cluster = '/output/%s/temp/%s' % (username_cluster, project_name)
		dir_cluster_ok = True
		
	while not input_dir_cluster_ok:
		print "Variable INPUT_DIR_CLUSTER not specified!"
		a = raw_input("INPUT_DIR_CLUSTER [/output/%s/temp/%s/input] = " % (username_cluster, project_name))
		if len(a) != 0:
			input_dir_cluster = copy.copy(a)
		else:
			input_dir_cluster = '/output/%s/temp/%s/input' % (username_cluster, project_name)
		input_dir_cluster_ok = True

	while not output_dir_cluster_ok:
		print "Variable OUTPUT_DIR_CLUSTER not specified!"
		a = raw_input("OUTPUT_DIR_CLUSTER [/output/%s/temp/%s/output] = " % (username_cluster, project_name))
		if len(a) != 0:
			output_dir_cluster = copy.copy(a)
		else:
			output_dir_cluster = '/output/%s/temp/%s/output' % (username_cluster, project_name)
		output_dir_cluster_ok = True

	while not scratch_dir_cluster_ok:
		print "Variable SCRATCH_DIR_CLUSTER not specified!"
		a = raw_input("SCRATCH_DIR_CLUSTER [/scratch/%s/%s] = " % (username_cluster, project_name))
		if len(a) != 0:
			scratch_dir_cluster = copy.copy(a)
		else:
			scratch_dir_cluster = '/scratch/%s/%s' % (username_cluster, project_name)
		scratch_dir_cluster_ok = True

	while not zindo_dir_cluster_ok:
		print "Variable ZINDO_DIR_CLUSTER not specified!"
		a = raw_input("ZINDO_DIR_CLUSTER [/home/output/David/Aijun/zindo-split/split] = ")
		if len(a) != 0:
			zindo_dir_cluster = copy.copy(a)
		else:
			zindo_dir_cluster = '/home/output/David/Aijun/zindo-split/split' % (username_cluster, project_name)
		zindo_dir_cluster_ok = True

	if quit:
		print "Aborting...\n"
		sys.exit(1)

	# ==========================
	#    Summary of variables
	# ==========================
	ClearScreen()

	summary = '===============================\n'
	summary += '      Summary of the input     \n'
	summary += '===============================\n\n'
	summary += 'PROJECT_NAME = %s\n' % (project_name)
	summary += 'INPUT_DIR = %s\n' % (input_dir)
	summary += 'FILE_TYPE = %s\n\n' % (file_type)
	summary += 'Neighbors parameters:\n'
	summary += '=====================\n'
	
	summary += 'MOLECULES_TO_ANALYZE = '
	for x in molecules_to_analyze:
		summary += ' %s ' % (x)
	summary += '\n'
		
	summary += 'PBC = '
	for x in pbc:
		summary += ' %s ' % (x)
	summary += '\n'
			
	summary += 'MOLECULES_FOR_J = '
	for x in molecules_for_J:
		summary += ' %s ' % (x)
	summary += '\n'
		
	summary += 'CUTOFF = %f \n\n' % (cutoff)
	summary += 'Cluster parameters:\n'
	summary += '===================\n'	
	summary += 'USERNAME_CLUSTER = %s \n' % (username_cluster)
	summary += 'DIR_CLUSTER = %s \n' % (dir_cluster)
	summary += 'INPUT_DIR_CLUSTER = %s \n' % (input_dir_cluster)
	summary += 'OUTPUT_DIR_CLUSTER = %s \n' % (output_dir_cluster)
	summary += 'SCRATCH_DIR_CLUSTER = %s \n' % (scratch_dir_cluster)
	summary += 'ZINDO_DIR_CLUSTER = %s \n' % (zindo_dir_cluster)
	summary += 'LOCATION_CLUSTER = %s \n' % (location_cluster)
	
	print summary
	
	# ==========================
	#   Check if informations
	# ==========================
	
	a = raw_input("Is it correct?[y,N] ")
	if a == "y" or a == "Y" or a == "yes":
		ClearScreen()
		try:
			finput = open('%s.sum' % (input_instructions), 'w')
			finput.write(summary)
			finput.close()
		except:
			if verb > 1:
				print "[WARNING] Could not write %s.sum for writing summary." % (input_instructions)
			pass
		if verb > 0:
			print "[YEAH] Informations are correct. In progress..."
	else:
		if verb > 0:
			print "[ERROR] You said that some informations are not correct. Aborting..."
		sys.exit(0)
		
	# ==========================
	#            PBC
	# ==========================
	pbc_number = ''
	if 'a' in pbc:
		pbc_number += '1 '
	else:
		pbc_number += '0 '

	if 'b' in pbc:
		pbc_number += '1 '
	else:
		pbc_number += '0 '
		
	if 'c' in pbc:
		pbc_number += '1'
	else:
		pbc_number += '0'
	pbc_number += '\n'
	
	# ==========================
	#     File type analysis
	# ==========================
	(file_type, ext) = FileTypeCheck(file_type)
	
	input_cluster = '%s%sinput%sMD' % (project_name, os.sep, os.sep)
	try:
		os.makedirs(input_cluster)
	except:
		pass
		if verb > 1:
			print "[WARNING] Could not create %s folder or folder already exists." % (input_cluster)
	
	try:
		shutil.copy("src%screate_input_zindo.cpp" % (os.sep), "%s" % (project_name))
	except:
		pass
		if verb > 1:
			print "[WARNING] Could not copy file src%create_input_zindo.cpp to %s.\n" % (os.sep, project_name)
	
	# ==========================
	#      Start working...
	# ==========================
	
	if verb > 2:
		print "[INFO] Work in progress..."
	for file in ls(input_dir, ext):
		# ==========================
		#     Input data analysis
		# ==========================
		if verb > 2:
			print "[INFO] Analyzing file %s" % (file)
		filename_base = os.path.splitext(file)[0]
		
		# Here where we create and fill the arrays.
		# Should be similar for other input files.
		if file_type == "tinker":
			file_add = "%s%s%s.add" % (input_dir, os.sep, filename_base)
			file_arc = "%s%s%s.arc" % (input_dir, os.sep, filename_base)
			
			# Getting informations about the size of the system and cell parameters
			(n_frame, n_mol, n_atom, a, b, c, alpha, beta, gamma) = input_reading_MD.Read_TINKER_add_File(file_add, verb)
			if verb > 2:
				print "[INFO] %s file read!" % (file_add)
			
			# Create MolecularSystem and SimulationBox
			box = molecular_system.SimulationBox(n_frame, a, b, c, alpha, beta, gamma)
			if verb > 2:
				print "[INFO] Generation of simulation box done!"
			qs = molecular_system.MolecularSystem(n_frame, n_mol, [n_atom])
			if verb > 2:
				print "[INFO] Generation of molecular system done!"

			# Reading Tinker file .arc
			input_reading_MD.Read_TINKER_arc_File(file_arc, qs, verb)
			if verb > 2:
				print "[INFO] %s file read!" % (file_arc)
			
			# Calculating quantities for the MolecularSystem and SimulationBox
			box.Parameters_For_Orthogonalization()
			if verb > 2:	
				print "[INFO] Parameters For Orthogonalization calculated!"
				print "[INFO] Calculation of Center of Masses. This may take some time..."
			qs.Center_of_Masses(list_manipulation.MoleculesList(molecules_to_analyze, qs.n_mol), verb)
			if verb > 2:	
				print "[INFO] Centers of Masses calculated!"
			
			# Making a link between qs and data
			data = qs
			cell = box
		
		# =============================
		# Include here other file types
		# =============================
		
		# ======================================
		# Analyzing the specified molecule lists
		# ======================================
		molecules_to_analyze_full = list_manipulation.MoleculesList(molecules_to_analyze, data.n_mol)
		molecules_for_J_full = list_manipulation.MoleculesList(molecules_for_J, data.n_mol)
		if verb > 2:
			print "[INFO] Creation of molecules lists done!"
		
		if(molecules_to_analyze != 'all' and molecules_for_J != 'all'):
			diff = list_manipulation.ListDifference(molecules_to_analyze_full, molecules_for_J_full)
			if len(diff) != 0:
				if verb > 0:
					print "[ERROR] Some of the molecules for J calculation are not in the\nlist of molecule considered for file %s:" % (file)
					print diff
					print "Aborting..."
				sys.exit(1)
				
			if int(molecules_to_analyze_full[len(molecules_to_analyze_full)-1]) > data.n_mol:
				if verb > 0:
					print "[ERROR] More molecules to analyze than molecules in the system for file %s." % (file)
					print "[ERROR] There are %d molecules in the system, and you asked for calculating molecule %d." % (data.n_mol, int(molecules_to_analyze_full[len(molecules_to_analyze_full)-1]))
					print "Aborting..."
				sys.exit(1)
			
		# ===================
		# Creating .xyz files
		# ===================
		
		if verb>2:
			print "[INFO] Writing output files. This may take some time..."
		cluster_file_creation.CreateXYZ(data, cell, molecules_to_analyze_full, input_cluster, filename_base, verb)
		cluster_file_creation.CreateCELL(data, cell, pbc_number, cutoff, input_cluster, filename_base)
		cluster_file_creation.CreateCM(data, molecules_to_analyze_full, molecules_for_J_full, input_cluster, filename_base)
		
	# =====================
	# Script files creation
	# =====================
	
	cluster_file_creation.ScriptFileCreation(dir_cluster, input_dir_cluster, output_dir_cluster, scratch_dir_cluster, zindo_dir_cluster, location_cluster, project_name)
	cluster_file_creation.ScriptFileCreationDirect(dir_cluster, input_dir_cluster, output_dir_cluster, zindo_dir_cluster, location_cluster, project_name)
	cluster_file_creation.ScriptFileCreationPBS(project_name, username_cluster, dir_cluster, location_cluster)
	cluster_file_creation.ScriptZINDOLaunch(dir_cluster, input_dir_cluster, output_dir_cluster, scratch_dir_cluster, project_name, location_cluster)

	try:
		import tarfile
		
		if verb > 2:
			print "[INFO] Generating compressed file. This may take some time..."
			
		tar = tarfile.open("%s.tar.gz" % project_name, "w:gz")
		tar.add("%s" % project_name)
		tar.close()
		
		if verb > 2:
			print "[INFO] Generation of compressed file done."
			
	except:
		if verb > 1:
			print "[WARNING] Could not create %s.tar.gz file.\n" % (project_name)	


	t2 = time.clock()
	if verb > 2:
		print "[INFO] The generation of the files took %f seconds.\n" % (t2-t1)
	
	sys.exit(0)
