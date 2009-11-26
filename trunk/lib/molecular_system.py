import math, numpy
import copy

#=====================================================================
#-------------------------Class Definition----------------------------
#=====================================================================

class SimulationBox(object):
	""" Class for the simulation box, containing the lattice parameters
	and temporary variables for fractional coordinates calculations """
	def __init__(self, n_frame, a, b, c, alpha, beta, gamma):
		self.n_frame = copy.copy(n_frame)
		
		self.a = copy.copy(a)
		self.b = copy.copy(b)
		self.c = copy.copy(c)
		
		self.alpha_deg = copy.copy(alpha)
		self.beta_deg = copy.copy(beta)
		self.gamma_deg = copy.copy(gamma)
		self.alpha_rad = numpy.zeros((n_frame), float)
		self.beta_rad = numpy.zeros((n_frame), float)
		self.gamma_rad = numpy.zeros((n_frame), float)
				
		for i in xrange(n_frame):
			self.alpha_rad[i] = math.radians(self.alpha_deg[i])
			self.beta_rad[i] = math.radians(self.beta_deg[i])
			self.gamma_rad[i] = math.radians(self.gamma_deg[i])
		
		self.temp_alpha_cos=numpy.zeros((n_frame), float)
		self.temp_beta_sin=numpy.zeros((n_frame), float)
		self.temp_beta_cos=numpy.zeros((n_frame), float)
		self.temp_gamma_sin=numpy.zeros((n_frame), float)
		self.temp_gamma_cos=numpy.zeros((n_frame), float)
		self.temp_beta_term=numpy.zeros((n_frame), float)
		self.temp_gamma_term=numpy.zeros((n_frame), float)

	def __str__(self):
		tmp = ''
		for i in xrange(self.n_frame):
			tmp += 'Frame %d\n' % (i)
			tmp += 'a=%f; b=%f; c=%f;\nalpha (deg)=%f; beta (deg)=%f; gamma (deg)=%f;\nalpha (rad)=%f; beta (rad)=%f; gamma (rad)=%f\n'\
			% (self.a[i], self.b[i], self.c[i], self.alpha_deg[i], self.beta_deg[i], self.gamma_deg[i], self.alpha_rad[i], self.beta_rad[i], self.gamma_rad[i])
		
		return tmp

#	def DegtoRad(self):
#		self.alpha_rad = math.radians(self.alpha_deg)
#		self.beta_rad = math.radians(self.beta_deg)
#		self.gamma_rad = math.radians(self.gamma_deg)
		
	def Parameters_For_Orthogonalization(self):
		for i in xrange(self.n_frame):
			if self.alpha_deg[i]==90.0 and self.beta_deg[i]==90.0 and self.gamma_deg[i]==90.0:
				self.temp_alpha_cos[i] = 0.0
				self.temp_beta_sin[i]= 1.0
				self.temp_beta_cos[i] = 0.0
				self.temp_gamma_sin[i] = 1.0
				self.temp_gamma_cos[i] = 0.0
				self.temp_beta_term[i] = 0.0
				self.temp_gamma_term[i] = 1.0

			elif self.alpha_deg[i]==90.0 and self.gamma_deg[i]==90.0:
				self.temp_alpha_cos[i] = 0.0
				self.temp_beta_sin[i] = math.sin(self.beta_rad[i])
				self.temp_beta_cos[i] = math.cos(self.beta_rad[i])
				self.temp_gamma_sin[i] = 1.0
				self.temp_gamma_cos[i] = 0.0
				self.temp_beta_term[i] = 0.0
				self.temp_gamma_term[i] = copy.copy(self.temp_beta_sin[i])

			else:
				self.temp_alpha_cos[i] = math.cos(self.alpha_rad[i])
				self.temp_beta_sin[i] = math.sin(self.beta_rad[i])
				self.temp_beta_cos[i] = math.cos(self.beta_rad[i])
				self.temp_gamma_sin[i] = math.sin(self.gamma_rad[i])
				self.temp_gamma_cos[i] = math.cos(self.gamma_rad[i])
				self.temp_beta_term[i] = (self.temp_alpha_cos[i]-self.temp_beta_cos[i]*self.temp_gamma_cos[i])/self.temp_gamma_sin[i]
				self.temp_gamma_term[i] = math.sqrt(1.0-math.pow(self.temp_beta_cos[i],2)-math.pow(self.temp_beta_term[i],2))
			
#		print self.temp_alpha_cos, self.temp_beta_sin, self.temp_beta_cos, self.temp_gamma_sin, self.temp_gamma_cos
#		print self.temp_beta_term, self.temp_gamma_term
		
class MolecularSystem(object):
	""" Class for atomic coordinates, containing also atomic symbol and masses
	and center of masses """
	def __init__(self, n_frame=1, n_mol=1, n_atom=[1]):
		self.n_frame = n_frame	# Number of frame
		self.n_mol = n_mol		# Number of molecules
		self.n_atom = numpy.zeros((n_mol), int)
		
		for ii in xrange(n_mol):
			if len(n_atom) != n_mol:
				self.n_atom[ii] = n_atom[0]
			else:
				self.n_atom[ii] = n_atom[ii]	# Number of atom per molecule 
		
		self.symbol = [[["A" for i in xrange(max(n_atom))] for j in xrange(n_mol)] for k in xrange(n_frame)]	# Atomic symbols
		self.x = numpy.zeros((n_frame, n_mol, max(n_atom)), float)	# x coordinates
		self.y = numpy.zeros((n_frame, n_mol, max(n_atom)), float)
		self.z = numpy.zeros((n_frame, n_mol, max(n_atom)), float)
		self.mol_number = numpy.zeros((n_frame, n_mol))	# Molecule label
		self.atomic_number = numpy.zeros((n_frame, n_mol, max(n_atom)), float)	# Atomic mass
		self.atomic_mass = numpy.zeros((n_frame, n_mol, max(n_atom)), float)	# Atomic mass
		self.atomic_valence = numpy.zeros((n_frame, n_mol, max(n_atom)), int)	# Atomic valence
		self.CM_x = numpy.zeros((n_frame, n_mol), float)	# Center of mass x coordinates
		self.CM_y = numpy.zeros((n_frame, n_mol), float)
		self.CM_z = numpy.zeros((n_frame, n_mol), float)
		self.n_electrons = numpy.zeros((n_frame, n_mol), int)
	
	def __str__(self):
		tmp = ""
		for i in xrange(self.n_frame):
			for ii in xrange(self.n_mol):
				for iii in xrange(self.n_atom[ii]):
					tmp += '%4s %10f %12f %12f %12f\n'\
					% (self.symbol[i][ii][iii], self.atomic_mass[i, ii, iii], self.x[i, ii, iii], self.y[i, ii, iii], self.z[i, ii, iii]) 
		return tmp
		
	def Center_of_Masses(self, molecules=[0], verb=2):
		if molecules[0] == 0:
			mol = xrange(self.n_mol)
		else:
			mol = molecules
			
		for i in xrange(self.n_frame):
			#for ii in xrange(self.n_mol):
			#for ii in molecules:
			for ii in mol:
				M_TOT = 0.0
				for iii in xrange(self.n_atom[ii]):
					self.n_electrons[i,ii] += self.atomic_valence[i, ii, iii]
					self.CM_x[i,ii] += self.atomic_mass[i, ii, iii] * self.x[i, ii, iii]
					self.CM_y[i,ii] += self.atomic_mass[i, ii, iii] * self.y[i, ii, iii]
					self.CM_z[i,ii] += self.atomic_mass[i, ii, iii] * self.z[i, ii, iii]
					M_TOT += self.atomic_mass[i, ii, iii]
				self.CM_x[i,ii] = self.CM_x[i,ii]/M_TOT
				self.CM_y[i,ii] = self.CM_y[i,ii]/M_TOT
				self.CM_z[i,ii] = self.CM_z[i,ii]/M_TOT
#				print M_TOT, self.CM_x[i,ii], self.CM_y[i,ii], self.CM_z[i,ii]
			if ( i%100 == 0 and verb>3):
				print "[INFO] CM of frame %d calculated!" % (i)

class Neighbors_System(object):
	""" Class with tables specifying if two molecules are neighbor_search.
	There is also the displacement vectors for PBC """
	def __init__(self, n_frame=1, n_mol_mol1=1, n_mol_mol2=1):
		self.n_frame = n_frame
		self.n_mol = n_mol_mol1
		self.neighbors = numpy.zeros((n_frame, n_mol_mol1, n_mol_mol2))
		self.displ_vec = numpy.zeros((n_frame, n_mol_mol1, n_mol_mol2, 3))
		self.dist_square = numpy.zeros((n_frame, n_mol_mol1, n_mol_mol2))

class Eigenvector_Matrix(object):
	""" Class for the eigenvectors of the matrix"""
	def __init__(self, n_modes=1):
		self.n_modes = n_modes		# Number of vibrational modes
		self.freq = numpy.zeros((self.n_modes), float)
		
		self.vec_matrix = numpy.matrix(numpy.zeros((self.n_modes, self.n_modes), float))	# Displacement vectors for each mode
		self.vec_matrix_I = numpy.matrix(numpy.zeros((self.n_modes, self.n_modes), float))	# Inverted matrix
		
		self.cart_coord_matrix = numpy.matrix(numpy.zeros((self.n_modes, 1), float))		# Cartesian coordinates matrix
		self.normal_coord_matrix = 0.0														# Cartesian coordinates matrix

	def __str__(self):
		tmp = ''
		for i in xrange(self.n_modes):
#		for i in xrange(10):
			tmp += 'Mode %d: %f cm-1\n' % (i+1, self.freq[i])
			for ii in xrange(0, self.n_modes, 3):
				tmp += '%6d %12.5f %12.5f %12.5f\n' % (ii/3 + 1, self.vec_matrix[i, ii], self.vec_matrix[i, ii+1], self.vec_matrix[i, ii+2])
			
		return tmp
		
	def getI(self):
		self.vec_matrix_I = self.vec_matrix.getI()
		
	def getRefCoord(self, mol):
		j = 0
		for i in xrange(mol.n_mol):
			for ii in xrange(mol.n_atom):
#				print i, ii, j
				self.cart_coord_matrix[j, 0] = mol.x[0,i,ii]
				self.cart_coord_matrix[j+1, 0] = mol.y[0,i,ii]
				self.cart_coord_matrix[j+2, 0] = mol.z[0,i,ii]
				j += 3
				
	def getNormalCoord(self):
		self.normal_coord_matrix = self.vec_matrix*self.cart_coord_matrix
