def isHMM():
    global msa
	for name,sequence in msa.items():
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