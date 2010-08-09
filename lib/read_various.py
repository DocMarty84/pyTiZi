import string
import os

def Read_Atomic_Numbers():
	""" Read atomic masses file """
	atomic_numbers = dict()
	
	try:
		finput = open("chemistry%sAtomic_Number" % (os.sep), 'r')
	except:
		print "Could not open chemistry%sAtomic_Number" % (os.sep)
		exit(1)
		
	for line in finput:
		words = line.split()
		words[0] = words[0].strip(string.punctuation + string.whitespace)
		atomic_numbers[words[0]] = float(words[1])
	return atomic_numbers

def Read_Atomic_Masses():
	""" Read atomic masses file """
	atomic_masses = dict()
	
	try:
		finput = open("chemistry%sAtomic_Masses" % (os.sep), 'r')
	except:
		print "Could not open chemistry%sAtomic_Masses" % (os.sep)
		exit(1)
		
	for line in finput:
		words = line.split()
		words[0] = words[0].strip(string.punctuation + string.whitespace)
		atomic_masses[words[0]] = float(words[1])
	return atomic_masses

def Read_Atomic_Valences():
	""" Read atomic valence file """
	atomic_valence = dict()
	
	try:
		finput = open("chemistry%sAtomic_Valence" % (os.sep), 'r')
	except:
		print "Could not open chemistry%sAtomic_Valence" % (os.sep)
		exit(1)
		
	for line in finput:
		words = line.split()
		words[0] = words[0].strip(string.punctuation + string.whitespace)
		atomic_valence[words[0]] = float(words[1])
	return atomic_valence
