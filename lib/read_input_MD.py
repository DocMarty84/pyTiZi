import string
import numpy as np
import sys
import read_various

def Read_PDB_File_First(name, verb=2):
	""" Calculate the number of frames/molecules/atoms in a PDB file,
	    and reads the box parameters. """
	try:
		finput = open(name, 'r')
	except:
		if verb > 0:
			print "[ERROR] Could not open %s" % (name)
		sys.exit(1)
	
	first_frame = True
	
	n_frame = 0
	n_mol = 0
	n_atom = []
	tmp = 0
	
	a_tmp = []
	b_tmp = []
	c_tmp = []
	alpha_tmp = []
	beta_tmp = []
	gamma_tmp = []
	
	while 1:
		line = finput.readline()
		words1 = line.split()
		
		if words1[0] == "END_OF_FILE":
			break
			
		elif words1[0] == "REMARK":
			pass
		
		elif words1[0] == "CRYST" or words1[0] == "CRYST1":
			a_tmp.append(float(words1[1]))
			b_tmp.append(float(words1[2]))
			c_tmp.append(float(words1[3]))
			alpha_tmp.append(float(words1[4]))
			beta_tmp.append(float(words1[5]))
			gamma_tmp.append(float(words1[6]))
			n_frame += 1
		
		elif (words1[0] == "ATOM" or words1[0] == "HETATM") and first_frame:
			tmp = 1
			n_mol = 1
			
			while True:
				line = finput.readline()
				words2 = line.split()
			
				if (words2[0] == "ATOM" or words2[0] == "HETATM") and int(words1[4]) == int(words2[4]):
					tmp += 1
					
				elif (words2[0] == "ATOM" or words2[0] == "HETATM") and int(words1[4]) != int(words2[4]):
					n_atom.append(tmp)
					words1 = words2[:]
					n_mol += 1
					tmp = 1
					
				elif (words2[0] == "TER" or words2[0] == "END"):
					n_atom.append(tmp)
					first_frame = False
					break

		else:
			pass
			
	finput.close()
	
	a = np.zeros((n_frame), float)
	b = np.zeros((n_frame), float)
	c = np.zeros((n_frame), float)
	alpha = np.zeros((n_frame), float)
	beta = np.zeros((n_frame), float)
	gamma = np.zeros((n_frame), float)
	
	for i in xrange(n_frame):
		a[i] = a_tmp[i]
		b[i] = b_tmp[i]
		c[i] = c_tmp[i]
		alpha[i] = alpha_tmp[i]
		beta[i] = beta_tmp[i]
		gamma[i] = gamma_tmp[i]
		
	return n_frame, n_mol, n_atom, a, b, c, alpha, beta, gamma

def Read_PDB_File_Second(name, mol, verb=2):
	""" Read the coordinates in a PDB file """
	atomic_numbers = read_various.Read_Atomic_Numbers()
	atomic_masses = read_various.Read_Atomic_Masses()
	atomic_valences = read_various.Read_Atomic_Valences()
	
	try:
		finput = open(name, 'r')
	except:
		if verb > 0:
			print "[ERROR] Could not open %s" % (name)
		sys.exit(1)
		
	i = -1
		
	while 1:		
		line = finput.readline()
		words = line.split()
		
		if words[0] == "END_OF_FILE":
			finput.close()
			break
			
		elif words[0] == "CRYST" or words[0] == "CRYST1":
			i += 1
			if verb > 3 and i%10 == 0:
				print "[INFO] Reading frame %d/%d" %(i, mol.n_frame)
			
		elif words[0] == "REMARK":
			pass
		
		elif words[0] == "ATOM" or words[0] == "HETATM":
			first_atom_line = True
			for ii in xrange(mol.n_mol):
				for iii in xrange(mol.n_atom[ii]):
					if not first_atom_line:
						line = finput.readline()
						words = line.split()
					mol.symbol[i][ii][iii] = words[2][0]
					mol.x[i, ii, iii] = float(words[5])
					mol.y[i, ii, iii] = float(words[6])
					mol.z[i, ii, iii] = float(words[7])
					mol.atomic_number[i, ii, iii] = atomic_numbers[mol.symbol[i][ii][iii]] 
					mol.atomic_mass[i, ii, iii] = atomic_masses[mol.symbol[i][ii][iii]] 
					mol.atomic_valence[i, ii, iii] = atomic_valences[mol.symbol[i][ii][iii]]
					
					first_atom_line = False
			
		else:
			pass
			
	finput.close()

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
