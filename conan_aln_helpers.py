import os
import conan_constants

def sortedSizeSequences():
    """
    Sorts sequences in the global MSA var from size (from biggest to smallest)
    """
    global msa
	new_seqs = []

	for seqname,sequence in msa.items():
		new_seq = sequence.replace('.','').replace('-','')
		new_seqs.append((seqname,new_seq))

	new_seqs = sorted(new_seqs,key=lambda x: len(x[1]),reverse=True)
	return new_seqs

def writeUnalignedFasta(fname="unaligned.fa", check=False):
    """
    From global MSA var, writes unaligned fasta file with name defined in fname variable (default: 'unaligned.fa')
    Input variable check verifies for existance of the fname variable named file before creating a new one
    """
    global msa
	unal_file = outputdir + "/" + fname
    if check and os.path.exists(unal_file):
            return unal_file
	fw = open(unal_file,'w')
	for seqname,sequence in msa.items():
		fw.write('>' + seqname + "\n")
		fw.write(sequence.replace('.','').replace('-','') + "\n")
	fw.close()
	return unal_file

def maxIdCDhit(check=False):
    """
    Returns new multiple sequence alignment filtered by the cd-hit program
    Cd-hit path should be defined in conan_constants.py

    Input variable check verifies for the 'cluster' file existance in the output directory before running cd-hit
    Before running Cd-hit, an unaligned fasta file is created (see above function)
    Actual sequences are retrieved from MSA by the 'cluster' file indexation
    """
    global msa
	unal_file = writeUnalignedFasta()
	out_file = outputdir + "/cluster"
	n = 5
	if maxid > 0.7:
		n = 5
	elif maxid > 0.6:
		n = 4
	elif maxid > 0.5:
		n = 3
	else:
		n = 2

	# os.system('cd-hit -i ' + unal_file + ' -o ' + out_file + ' -c ' + str(maxid) + ' -n ' + str(n) + ' -M 3000 -T 2')

	if (not check) or (check and os.path.exists(out_file) == False):
        os.system(conan_constants.CD_HIT_PATH + ' -i ' + unal_file + ' -o ' + out_file + ' -c ' + str(maxid) + ' -n ' + str(n) + ' -M 3000 -T 2')

	new_msa = {}
	fr = open(out_file)

	for line in fr:
		line = line.strip()
		if len(line) > 1:
			if line[0] == '>':
				seqname = line[1:]
                #not retrieved from file but from msa
				sequence = msa[seqname]
				new_msa[seqname] = sequence
	fr.close()
	return new_msa

def aa2id(aa):
    """
    Converts Aminoacids ACDEFGHIKLMNPQRSTVWY to indexa 0-19 and gaps -. to 20
    All else becomes index 21
    """
	aa = aa.upper()
	if aa == 'A':
		return 0
	elif aa == 'C':
		return 1
	elif aa == 'D':
		return 2
	elif aa == 'E':
		return 3
	elif aa == 'F':
		return 4
	elif aa == 'G':
		return 5
	elif aa == 'H':
		return 6
	elif aa == 'I':
		return 7
	elif aa == 'K':
		return 8
	elif aa == 'L':
		return 9
	elif aa == 'M':
		return 10
	elif aa == 'N':
		return 11
	elif aa == 'P':
		return 12
	elif aa == 'Q':
		return 13
	elif aa == 'R':
		return 14
	elif aa == 'S':
		return 15
	elif aa == 'T':
		return 16
	elif aa == 'V':
		return 17
	elif aa == 'W':
		return 18
	elif aa == 'Y':
		return 19
	elif aa == '-':
		return 20
	elif aa == '.':
		return 20
	else:
		return 21

def readFasta(inputfile):
    """
    Read Fasta formatted MSA file, outputs dict with header : sequence
    """
	msa = {}
	fr = open(inputfile,'r')
	header = None
	sequence = ""

	for line in fr:
		line = line.strip()
		if line.startswith('>'):
			if header is not None:
				msa[header] = sequence
				header = line[1:]
				sequence = ""
			else:
				header = line[1:]
		else:
			sequence += line
	if header is not None:
		msa[header] = sequence

	fr.close()
	return msa

def readStockholm(inputfile):
    """
    Read Pfam/Stockholm formatted MSA file, outputs dict with header : sequence
    """
	msa = {}
	fr = open(inputfile,'r')

	for line in fr:
		if line[0] != "#":
			line = line.strip()
			temp = line.split()
			if len(temp) > 1:
				seqname = temp[0]
				sequence = temp[1]
				msa[seqname] = sequence

	fr.close()
	return msa