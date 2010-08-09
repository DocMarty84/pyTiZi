import string

def ListDifference(a, b):
	""" show whats in list b which isn't in list a """
	return list(set(b).difference(set(a))) 

def Make1DList(alist):
	return [ int(item) for sublist in alist for item in sublist ]

def RemoveDuplicate(alist): # Fastest order preserving
	set = {}
	return [set.setdefault(e,e) for e in alist if e not in set]

def MoleculesList(molecules, n_mol):
	if molecules[0] == 'all':
		molecules_full = range(n_mol)
		
	else:
		a = []
		
		for x in molecules:
			n = x.split("-")
			
			# If molecule alone
			if len(n) == 1:
				a.append(n)
			
			# If something like 3-7 is specified
			elif len(n) == 2:
				if int(n[0]) < int(n[1]):
					t_min = int(n[0])
					t_max = int(n[1]) + 1
					m = range(t_min, t_max, 1)
				elif int(n[0]) > int(n[1]):
					t_min = int(n[1])
					t_max = int(n[0]) + 1
					m = range(t_min, t_max, 1)
				else:
					m = int(n[0])
				a.append(m)
				
			else:
				print "Check the variable MOLECULES_TO_ANALYZE, it seems that there is a problem..."
		
		# Make a 1D list
		b = Make1DList(a)

		# Remove duplicates and sort
		molecules_full = RemoveDuplicate(b)
		molecules_full.sort()

	return molecules_full
