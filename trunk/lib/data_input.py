import os
import re

class InstructionsData(object):
	""" Class containing all intruction for the simulation """
	def __init__(self):
		self.project_name = ''
		self.input_dir = ''
		self.file_type = ''
		self.input_dir_local = ''
		self.molecules_to_analyze = ''
		self.pbc = ''
		self.molecules_for_J = ''
		self.cutoff = ''
		self.username_cluster = ''
		self.dir_cluster = ''
		self.input_dir_cluster = ''
		self.output_dir_cluster = ''
		self.results_dir_cluster = ''
		self.scratch_dir_cluster = ''
		self.location_cluster = ''
		self.zindo_dir_cluster = ''
		self.sign = ''
		self.coeff_H_lign = 0
		self.coeff_H_row = 0
		self.coeff_L_lign = 0
		self.coeff_L_row = 0
		
		self.project_name_ok = False
		self.input_dir_ok = False
		self.file_type_ok = False
		self.input_dir_local_ok = False
		self.molecules_to_analyze_ok = False
		self.pbc_ok = False
		self.molecules_for_J_ok = False
		self.cutoff_ok = False
		self.username_cluster_ok = False
		self.dir_cluster_ok = False
		self.input_dir_cluster_ok = False
		self.output_dir_cluster_ok = False
		self.results_dir_cluster_ok = False
		self.scratch_dir_cluster_ok = False
		self.location_cluster_ok = False
		self.zindo_dir_cluster_ok = False
		self.sign_ok = False
		
		self.input_cluster = ''
		self.pbc_number = ''
		self.ext = ''
		self.molecules_to_analyze_full = ''
		self.molecules_for_J_full = ''

	def __str__(self):
		summary = '===============================\n'
		summary += '      Summary of the input     \n'
		summary += '===============================\n\n'
		summary += 'PROJECT_NAME = %s\n' % (self.project_name)
		summary += 'INPUT_DIR = %s\n' % (self.input_dir)
		summary += 'FILE_TYPE = %s\n\n' % (self.file_type)
		summary += 'Neighbors parameters:\n'
		summary += '=====================\n'
		
		summary += 'MOLECULES_TO_ANALYZE = '
		for x in self.molecules_to_analyze:
			summary += ' %s ' % (x)
		summary += '\n'
			
		summary += 'PBC = '
		for x in self.pbc:
			summary += ' %s ' % (x)
		summary += '\n'
				
		summary += 'MOLECULES_FOR_J = '
		for x in self.molecules_for_J:
			summary += ' %s ' % (x)
		summary += '\n'
		
		summary += 'CALCULATE_SIGN = %d, %d %d %d %d \n' % (self.sign, self.coeff_H_lign, self.coeff_H_row, self.coeff_L_lign, self.coeff_L_row)

			
		summary += 'CUTOFF = %f \n\n' % (self.cutoff)
		summary += 'Cluster parameters:\n'
		summary += '===================\n'	
		summary += 'USERNAME_CLUSTER = %s \n' % (self.username_cluster)
		summary += 'DIR_CLUSTER = %s \n' % (self.dir_cluster)
		summary += 'SCRATCH_DIR_CLUSTER = %s \n' % (self.scratch_dir_cluster)
		summary += 'ZINDO_DIR_CLUSTER = %s \n' % (self.zindo_dir_cluster)
		summary += 'LOCATION_CLUSTER = %s \n' % (self.location_cluster)
		
		return summary


def ls(path, ext):
	return [fichier for fichier in os.listdir(path)
								if os.path.splitext(fichier)[1] == ext]
								
def grep(string,list):
	expr = re.compile(string)
	return filter(expr.search,list)
		
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
