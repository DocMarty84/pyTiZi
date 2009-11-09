#!/usr/bin/env python

import string
import os, shutil, sys
import time
import copy

sys.path.append("lib")
import read_input_MD, molecular_system, list_manipulation, write_cluster_files
from data_input import *
			
if __name__ == '__main__':
	t1 = time.clock()
	
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

	if len(sys.argv) > 1:
		tmp = sys.argv[1:]
		input_instructions = ' '.join(x for x in tmp)
		
	else:
		try:
			import Tkinter
			import tkFileDialog
			input_instructions = tkFileDialog.askopenfilename(filetypes = [("All", "*")], title='Please select an instruction file')
		except ImportError:
			if verb > 1:
				print "[WARNING] Tkinter or tkFileDialog could not be imported. Instruction file should be named \"test\"."
				input_instructions = "test"

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
	#   Instruction data class
	# ==========================
	
	project = InstructionsData()
	
	# ==========================
	# Reading instructions file
	# ==========================
	try:
		finstructions = open(input_instructions, 'r')
	except:
		if verb > 0:
			print "[ERROR] Could not find %s file. Aborting.\n" % (input_instructions)
		sys.exit(1)
	
	for line in finstructions:
		words = line.split()
		if len(words) != 0:
			words[0] = words[0].strip(string.punctuation + string.whitespace)
			first_word = string.upper(words[0])
			
			if first_word=="PROJECT_NAME":
				project.project_name = words[1]
				project.project_name_ok = True
			
			elif first_word=="INPUT_DIR":
				tmp = words[1:]
				project.input_dir = ' '.join(x for x in tmp)
				project.input_dir_ok = True
				
			elif first_word=="FILE_TYPE":
				project.file_type = string.lower(words[1])
				project.file_type_ok = True
			
			elif first_word=="MOLECULES_TO_ANALYZE":
	#			cut = re.compile('[;, ]')
	#			molecules_to_analyze = [ cut.split(x) for x in words[1:] ]
				a = [ x.split(",") for x in words[1:] ]
				project.molecules_to_analyze = [ string.lower(item) for sublist in a for item in sublist ]
				project.molecules_to_analyze_ok = True
			
			elif first_word=="PBC":
				a = [ x.split(",") for x in words[1:] ]
				project.pbc = [ string.lower(item) for sublist in a for item in sublist ]
				project.pbc_ok = True
			
			elif first_word=="MOLECULES_FOR_J":
				a = [ x.split(",") for x in words[1:] ]
				project.molecules_for_J = [ string.lower(item) for sublist in a for item in sublist ]
				project.molecules_for_J_ok = True
			
			elif first_word=="CUTOFF":
				project.cutoff = float(words[1])
				project.cutoff_ok = True
			
			elif first_word=="USERNAME_CLUSTER":
				project.username_cluster = words[1]
				project.username_cluster_ok = True

			elif first_word=="DIR_CLUSTER":
				project.dir_cluster = words[1]
				project.dir_cluster_ok = True

			elif first_word=="INPUT_DIR_CLUSTER":
				project.input_dir_cluster = words[1]
				project.input_dir_cluster_ok = True
				
			elif first_word=="OUTPUT_DIR_CLUSTER":
				project.output_dir_cluster = words[1]
				project.output_dir_cluster_ok = True
			
			elif first_word=="SCRATCH_DIR_CLUSTER":
				project.scratch_dir_cluster = words[1]
				project.scratch_dir_cluster_ok = True
				
			elif first_word=="LOCATION_CLUSTER":
				project.location_cluster = words[1]
				project.location_cluster_ok = True
			
			elif first_word=="ZINDO_DIR_CLUSTER":
				project.zindo_dir_cluster = words[1]
				project.zindo_dir_cluster_ok = True
				
			elif first_word=="END":
				break

	finstructions.close()

	quit = False
	
	# ==========================
	#   Variables obligatoires
	# ==========================
	
	while not project.project_name_ok:
		print "Variable PROJECT_NAME not specified!"
		a = raw_input("PROJECT_NAME = ")
		while len(a) == 0:
			a = raw_input("PROJECT_NAME = ")
			
		project.project_name = copy.copy(a)
		project.project_name_ok = True

	while not project.file_type_ok:
		print "Variable FILE_TYPE not specified!"
		a = raw_input("FILE_TYPE = ")
		while len(a) == 0:
			a = raw_input("FILE_TYPE = ")
			
		project.file_type = string.lower(a)
		project.file_type_ok = True
	
	while not project.molecules_to_analyze_ok:
		print "Variable MOLECULES_TO_ANALYZE not specified!"
		a = raw_input("MOLECULES_TO_ANALYZE = ")
		while len(a) == 0:
			a = raw_input("MOLECULES_TO_ANALYZE = ")
		
		a = [ x.split(",") for x in a ]
		project.molecules_to_analyze = [ string.lower(item) for sublist in a for item in sublist ]
		project.molecules_to_analyze_ok = True
		
	while not project.pbc_ok:
		print "Variable PBC not specified!"
		a = raw_input("PBC = ")
		while len(a) == 0:
			a = raw_input("PBC = ")

		a = [ x.split(",") for x in a ]
		project.pbc = [ string.lower(item) for sublist in a for item in sublist ]
		project.pbc_ok = True

	while not project.molecules_for_J_ok:
		print "Variable MOLECULES_FOR_J not specified!"
		a = raw_input("MOLECULES_FOR_J = ")
		while len(a) == 0:
			a = raw_input("MOLECULES_FOR_J = ")

		a = [ x.split(",") for x in a ]
		project.molecules_for_J = [ string.lower(item) for sublist in a for item in sublist ]
		project.molecules_for_J_ok = True

	while not project.cutoff_ok:
		print "Variable CUTOFF not specified!"
		a = raw_input("CUTOFF = ")
		while len(a) == 0:
			a = raw_input("CUTOFF = ")
			
		project.cutoff = float(a)
		project.cutoff_ok = True

	while not project.username_cluster_ok:
		print "Variable USERNAME_CLUSTER not specified!"
		a = raw_input("USERNAME_CLUSTER = ")
		while len(a) == 0:
			a = raw_input("USERNAME_CLUSTER = ")
			
		project.username_cluster = copy.copy(a)
		project.username_cluster_ok = True
		
	while not project.location_cluster_ok:
		print "Variable LOCATION_CLUSTER not specified!"
		a = raw_input("LOCATION_CLUSTER = ")
		while len(a) == 0:
			a = raw_input("LOCATION_CLUSTER = ")
			
		project.location_cluster = copy.copy(a)
		project.location_cluster_ok = True

	# ==========================
	# Variables non obligatoires
	# ==========================

	while not project.input_dir_ok:
		print "Variable INPUT_DIR not specified!"
		a = raw_input("INPUT_DIR [input_data] = ")
		if len(a) != 0:
			project.input_dir = copy.copy(a)
		else:
			project.input_dir = "input_data"
		project.input_dir_ok = True

	while not project.dir_cluster_ok:
		print "Variable DIR_CLUSTER not specified!"
		a = raw_input("DIR_CLUSTER [/output/%s/temp/%s] = " % (project.username_cluster, project.project_name))
		if len(a) != 0:
			project.dir_cluster = copy.copy(a)
		else:
			project.dir_cluster = '/output/%s/temp/%s' % (project.username_cluster, project.project_name)
		project.dir_cluster_ok = True
		
	while not project.input_dir_cluster_ok:
		print "Variable INPUT_DIR_CLUSTER not specified!"
		a = raw_input("INPUT_DIR_CLUSTER [/output/%s/temp/%s/input] = " % (project.username_cluster, project.project_name))
		if len(a) != 0:
			project.input_dir_cluster = copy.copy(a)
		else:
			project.input_dir_cluster = '/output/%s/temp/%s/input' % (project.username_cluster, project.project_name)
		project.input_dir_cluster_ok = True

	while not project.output_dir_cluster_ok:
		print "Variable OUTPUT_DIR_CLUSTER not specified!"
		a = raw_input("OUTPUT_DIR_CLUSTER [/output/%s/temp/%s/output] = " % (project.username_cluster, project.project_name))
		if len(a) != 0:
			project.output_dir_cluster = copy.copy(a)
		else:
			project.output_dir_cluster = '/output/%s/temp/%s/output' % (project.username_cluster, project.project_name)
		project.output_dir_cluster_ok = True

	while not project.scratch_dir_cluster_ok:
		print "Variable SCRATCH_DIR_CLUSTER not specified!"
		a = raw_input("SCRATCH_DIR_CLUSTER [/scratch/%s/%s] = " % (project.username_cluster, project.project_name))
		if len(a) != 0:
			project.scratch_dir_cluster = copy.copy(a)
		else:
			project.scratch_dir_cluster = '/scratch/%s/%s' % (project.username_cluster, project.project_name)
		project.scratch_dir_cluster_ok = True

	while not project.zindo_dir_cluster_ok:
		print "Variable ZINDO_DIR_CLUSTER not specified!"
		a = raw_input("ZINDO_DIR_CLUSTER [/home/output/David/Aijun/zindo-split/split] = ")
		if len(a) != 0:
			project.zindo_dir_cluster = copy.copy(a)
		else:
			project.zindo_dir_cluster = '/home/output/David/Aijun/zindo-split/split' % (project.username_cluster, project.project_name)
		project.zindo_dir_cluster_ok = True

	if quit:
		print "Aborting...\n"
		sys.exit(1)

	# ==========================
	#    Summary of variables
	# ==========================
	ClearScreen()
	print project
	
	# ==========================
	#   Check if informations
	# ==========================
	
	a = raw_input("Is it correct?[y,N] ")
	if a == "y" or a == "Y" or a == "yes":
		ClearScreen()
		try:
			finput = open('%s.sum' % (input_instructions), 'w')
			tmp = project.__str__()
			finput.write(tmp)
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
	project.pbc_number = ''
	if 'a' in project.pbc:
		project.pbc_number += '1 '
	else:
		project.pbc_number += '0 '

	if 'b' in project.pbc:
		project.pbc_number += '1 '
	else:
		project.pbc_number += '0 '
		
	if 'c' in project.pbc:
		project.pbc_number += '1'
	else:
		project.pbc_number += '0'
	project.pbc_number += '\n'
	
	# ==========================
	#     File type analysis
	# ==========================
	(project.file_type, project.ext) = FileTypeCheck(project.file_type)
	
	project.input_cluster = 'project%s%s%sinput%sMD' % (os.sep, project.project_name, os.sep, os.sep)
	try:
		os.makedirs(project.input_cluster)
	except:
		pass
		if verb > 1:
			print "[WARNING] Could not create %s folder or folder already exists." % (project.input_cluster)
	
	try:
		shutil.copy("src%screate_input_zindo.cpp" % (os.sep), "project%s%s" % (os.sep, project.project_name))
	except:
		pass
		if verb > 1:
			print "[WARNING] Could not copy file src%create_input_zindo.cpp to %s.\n" % (os.sep, project.project_name)
	
	# ==========================
	#      Start working...
	# ==========================
	
	if verb > 2:
		print "[INFO] Work in progress..."
	for file in ls(project.input_dir, project.ext):
		# ==========================
		#     Input data analysis
		# ==========================
		if verb > 2:
			print "[INFO] Analyzing file %s" % (file)
		filename_base = os.path.splitext(file)[0]
		
		# Here where we create and fill the arrays.
		# Should be similar for other input files.
		if project.file_type == "tinker":
			file_add = "%s%s%s.add" % (project.input_dir, os.sep, filename_base)
			file_arc = "%s%s%s.arc" % (project.input_dir, os.sep, filename_base)
			
			# Getting informations about the size of the system and cell parameters
			(n_frame, n_mol, n_atom, a, b, c, alpha, beta, gamma) = read_input_MD.Read_TINKER_add_File(file_add, verb)
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
			read_input_MD.Read_TINKER_arc_File(file_arc, qs, verb)
			if verb > 2:
				print "[INFO] %s file read!" % (file_arc)
			
			# Calculating quantities for the MolecularSystem and SimulationBox
			box.Parameters_For_Orthogonalization()
			if verb > 2:	
				print "[INFO] Parameters For Orthogonalization calculated!"
				print "[INFO] Calculation of Center of Masses. This may take some time..."
			qs.Center_of_Masses(list_manipulation.MoleculesList(project.molecules_to_analyze, qs.n_mol), verb)
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
		project.molecules_to_analyze_full = list_manipulation.MoleculesList(project.molecules_to_analyze, data.n_mol)
		project.molecules_for_J_full = list_manipulation.MoleculesList(project.molecules_for_J, data.n_mol)
		if verb > 2:
			print "[INFO] Creation of molecules lists done!"
		
		if(project.molecules_to_analyze != 'all' and project.molecules_for_J != 'all'):
			diff = list_manipulation.ListDifference(project.molecules_to_analyze_full, project.molecules_for_J_full)
			if len(diff) != 0:
				if verb > 0:
					print "[ERROR] Some of the molecules for J calculation are not in the\nlist of molecule considered for file %s:" % (file)
					print diff
					print "Aborting..."
				sys.exit(1)
				
			if int(project.molecules_to_analyze_full[len(project.molecules_to_analyze_full)-1]) > data.n_mol:
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
		write_cluster_files.CreateXYZ(data, cell, project, filename_base, verb)
		write_cluster_files.CreateCELL(data, cell, project, filename_base)
		write_cluster_files.CreateCM(data, project, filename_base)
		
	# =====================
	# Script files creation
	# =====================
	
	write_cluster_files.ScriptFileCreation(project)
	write_cluster_files.ScriptFileCreationDirect(project)
	write_cluster_files.ScriptFileCreationPBS(project)
	write_cluster_files.ScriptZINDOLaunch(project)

	try:
		import tarfile
		
		if verb > 2:
			print "[INFO] Generating compressed file. This may take some time..."
		
		os.chdir("project")	
		tar = tarfile.open("%s.tar.gz" % (project.project_name), "w:gz")
		tar.add("%s" % (project.project_name))
		tar.close()
		
		if verb > 2:
			print "[INFO] Generation of compressed file done."
			
	except:
		if verb > 1:
			print "[WARNING] Could not create project%s%s.tar.gz file.\n" % (os.sep, project.project_name)	


	t2 = time.clock()
	if verb > 2:
		print "[INFO] The generation of the files took %f seconds.\n" % (t2-t1)
	
	sys.exit(0)
