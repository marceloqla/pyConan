import shared_globals

def isHMM():
	for name,sequence in shared_globals.msa.items():
		if '.' in sequence:
			return True
	return False	

def isAA(aa,case):
	if case:
		if aa != "." and aa != "-" and aa.isupper():
			return True
		else:
			return False
	else:
		if aa != "." and aa != "-":
			return True
		else:
			return False