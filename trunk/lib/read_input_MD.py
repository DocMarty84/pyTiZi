import string
import numpy as np
import sys
import read_various

def Calculate_Number_Atom_Mol_PDB(name):
	""" Calculate the number of different molecules in a pdb file """
	try:
		finput = open(name, 'r')
	except:
		print "Could not open %s" % (name)
		sys.exit(1)
	
	n_mol = 0
	n_atom = []
	tmp = 0
	
	while 1:
		line = finput.readline()
		words1 = line.split()
		if words1[0] == "ATOM":
			line = finput.readline()
			words2 = line.split()
			
			while words2[0] != "TER" and int(words1[4]) == int(words2[4]):
				tmp += 1
				line = finput.readline()
				words2 = line.split()
				
			tmp += 1	
			n_mol += 1
			n_atom.append(tmp)
			tmp = 1
			
			if words2[0] == "TER":
				break
	
	return n_mol, n_atom

def Read_PDB_Interface(name, mol1, mol2):
	""" Read a PDB file (suppose 2 kind of molecules) """
	atomic_masses = Read_Atomic_Masses()
	
	try:
		finput = open(name, 'r')
	except:
		print "Could not open %s" % (name)
		sys. exit(1)
		
	while 1:		
		line = finput.readline()
		words = line.split()
		first_word = words[0].strip(string.punctuation + string.whitespace)
		if first_word == "REMARK":
			pass

		elif first_word == "CRYST":
			box.a = float(words[1])
			box.b = float(words[2])
			box.c = float(words[3])
			box.alpha_deg = float(words[4])
			box.beta_deg = float(words[5])
			box.gamma_deg = float(words[6])
			box.DegtoRad()

		elif first_word == "ATOM":
			first_atom_line = True
			for i in xrange(mol1.n_frame):
				for ii in xrange(mol1.n_mol):
					for iii in xrange(mol1.n_atom[ii]):
						if not first_atom_line:
							line = finput.readline()
							words = line.split()
						mol1.symbol[i][ii][iii] = words[2].strip(string.punctuation + string.whitespace)
						mol1.x[i, ii, iii] = float(words[5])
						mol1.y[i, ii, iii] = float(words[6])
						mol1.z[i, ii, iii] = float(words[7])
						mol1.atomic_mass[i, ii, iii] = atomic_masses[mol1.symbol[i][ii][iii]]
						first_atom_line = False
#						print mol1.symbol[i][ii][iii], mol1.atomic_mass[i, ii, iii], mol1.x[i, ii, iii], mol1.y[i, ii, iii], mol1.z[i, ii, iii]
					mol1.mol_number[i, ii] = int(words[4])
					
			for j in xrange(mol2.n_frame):
				for jj in xrange(mol2.n_mol):
					for jjj in xrange(mol2.n_atom[jj]):
						line = finput.readline()
						words = line.split()
						mol2.symbol[j][jj][jjj] = words[2].strip(string.punctuation + string.whitespace)
						mol2.x[j, jj, jjj] = float(words[5])
						mol2.y[j, jj, jjj] = float(words[6])
						mol2.z[j, jj, jjj] = float(words[7])
						mol2.atomic_mass[j, jj, jjj] = atomic_masses[mol2.symbol[j][jj][jjj]]
#						print mol2.symbol[j][jj][jjj], mol2.atomic_mass[j, jj, jjj], mol2.x[j, jj, jjj], mol2.y[j, jj, jjj], mol2.z[j, jj, jjj]
					mol2.mol_number[j, jj] = int(words[4])
		elif first_word == "TER":
			break
			
	finput.close()
	print "Reading of input file done!"

def Read_TINKER_add_File(name, verb=2):
	try:
		finput = open(name, 'r')
	except:
		if verb>0:
			print "[ERROR] Could not open %s. Aborting..." % (name)
		sys.exit(1)
	
	for line in finput:
		words = line.split()
		if len(words) != 0:
			words[0] = words[0].strip(string.punctuation + string.whitespace)
			first_word = string.lower(words[0])
			
			if first_word == "simple":
				pass
			
			elif first_word == "n_frame":
				n_frame = int(words[1])
				
				a = np.zeros((n_frame), float)
				b = np.zeros((n_frame), float)
				c = np.zeros((n_frame), float)
				alpha = np.zeros((n_frame), float)
				beta = np.zeros((n_frame), float)
				gamma = np.zeros((n_frame), float)
				
			elif first_word == "n_mol":
				n_mol = int(words[1])
				
			elif first_word == "n_atom":
				n_atom = int(words[1])
			
			elif first_word == "frame":
				frame = int(words[1])
				if frame+1 > n_frame:
					if verb > 0:
						print "[ERROR] While reading %s: at least one frame \
						       number is higher than the total number of \
						       frames. Aborting..." % (name)
						sys.exit(1)
						
				for i in xrange(2, 14, 1):
					if words[i] == "a":
						a[frame] = float(words[i+1])
						
					if words[i] == "b":
						b[frame] = float(words[i+1])
											
					if words[i] == "c":
						c[frame] = float(words[i+1])
						
					if words[i] == "alpha":
						alpha[frame] = float(words[i+1])
						
					if words[i] == "beta":
						beta[frame] = float(words[i+1])
						
					if words[i] == "gamma":
						gamma[frame] = float(words[i+1])
						
	finput.close()
	
	return (n_frame, n_mol, n_atom, a, b, c, alpha, beta, gamma)
	
def Read_TINKER_arc_File(name, mol, verb=2):
	atomic_numbers = read_various.Read_Atomic_Numbers()
	atomic_masses = read_various.Read_Atomic_Masses()
	atomic_valences = read_various.Read_Atomic_Valences()
	try:
		finput = open(name, 'r')
	except:
		if verb > 0:
			print "[ERROR] Could not open %s. Aborting..." % (name)
		sys.exit(1)
	
	for i in xrange(mol.n_frame):
		line = finput.readline()
		words = line.split()
		while len(words) == 0:
			line = finput.readline()
			words = line.split()
		
		for ii in xrange(mol.n_mol):
			for iii in xrange(mol.n_atom[ii]):
				line = finput.readline()
				words = line.split()
				mol.symbol[i][ii][iii] = words[1]
				mol.x[i, ii, iii] = float(words[2])
				mol.y[i, ii, iii] = float(words[3])
				mol.z[i, ii, iii] = float(words[4])
				mol.atomic_number[i, ii, iii] = atomic_numbers[mol.symbol[i][ii][iii]] 
				mol.atomic_mass[i, ii, iii] = atomic_masses[mol.symbol[i][ii][iii]] 
				mol.atomic_valence[i, ii, iii] = atomic_valences[mol.symbol[i][ii][iii]] 
		if ( i%100 == 0 and verb > 3):
			print "[INFO] Frame %d read!" % (i)

def Read_TINKER_vec_File(name, mol, data):
	""" Read Tinker Vector file """
	try:
		finput = open(name, 'r')
	except:
		print "Could not open %s" % (name)
		exit(1)

	i = 0
	
	line = finput.readline()

	while len(line)!=0:
#			print len(line)
		words = line.split()
		
		while (len(words)==0):
			line = finput.readline()
			words = line.split()
			
		if (len(words)!=0 and (words[0]=="Vibrational" and words[1]=="Normal" and words[2]=="Mode")):
			try:
				data.freq[i] = float(words[6])
			except:
				data.freq[i] = 0.0
				
			line = finput.readline()
			words = line.split()
			
			while (len(words)==0):
				line = finput.readline()
				words = line.split()
				
			if (len(words)!=0 and (words[0]=="Atom" and words[1]=="Delta" and words[2]=="X")):
				line = finput.readline()
				words = line.split()
				
				while (len(words)==0):
					line = finput.readline()
					words = line.split()
				
				data.vec_matrix[i,0] = float(words[1])
				data.vec_matrix[i,1] = float(words[2])
				data.vec_matrix[i,2] = float(words[3])
				for ii in xrange(1, mol.n_mol * mol.n_atom[0], 1):
					line = finput.readline()
					words = line.split()
					
					data.vec_matrix[i,ii*3] = float(words[1])
					data.vec_matrix[i,(ii*3) + 1] = float(words[2])
					data.vec_matrix[i,(ii*3) + 2] = float(words[3])
					
			i += 1
		line = finput.readline()
		
	finput.close()
		
#	print "Reading of input file done!"
