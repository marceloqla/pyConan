from scipy.stats import hypergeom
from scipy.special import logsumexp
import conan_constants

def pearsonEdge(res1,res2):
    """
    WARN: Currently unused
    Calculate Pearson edge correlation for two residue Node objects (res1, res2)
    """
    global msa

	aa1_list = conan_constants.setsAA[res1.id[0]]
	i1 = res1.position[0]-1
	aa2_list = conan_constants.setsAA[res2.id[0]]
	i2 = res2.position[0]-1
	sx = 0.0
	sy = 0.0
	sxx = 0.0
	syy = 0.0
	sxy = 0.0
	n = float(len(msa))

	for sequence in msa.values():
		x = 0.0
		y = 0.0

		if sequence[i1].upper() in aa1_list:
			x=1.0
		if sequence[i2].upper() in aa2_list:
			y=1.0
		sx += x
		sy += y
		sxx += x * x
		syy += y * y
		sxy += x * y

	#print(str(aa1_list) + " " + str(aa2_list) + " " + str(sx) + " " + str(sy) + " " + str(sxx) + " " + str(syy) + " " + str(sxy) + " " + str(n))
	cov = sxy / n - sx * sy / n / n
	sigmax = math.sqrt(sxx / n -  sx * sx / n / n)
	sigmay = math.sqrt(syy / n -  sy * sy / n / n)
	return cov / sigmax / sigmay

def borgatti(res1,res2,type):
    """
    WARN: Currently unused
    Calculate Borgatti edge correlation for two residue Node objects (res1, res2)
    """
    global msa

	aa1_list = conan_constants.setsAA[res1.id[0]]
	i1 = res1.position[0]-1
	aa2_list = conan_constants.setsAA[res2.id[0]]
	i2 = res2.position[0]-1
	a = 0.0
	b = 0.0
	c = 0.0
	d = 0.0

	for sequence in msa.values():
		bool1 = sequence[i1].upper() in aa1_list
		bool2 = sequence[i2].upper() in aa2_list

		if bool1 and bool2:
			a+=1.0
		elif bool1:
			b += 1.0
		elif bool2:
			c += 1.0
		else:
			d += 1.0

	if type == 0:#BN
		if a+b > a+c:
			return a/(a+c)
		else:
			return a/(a+b)
	elif type == 1:#JC
		return a/(a+b+c)
	elif type == 2:#Bonacich
		if a*d == b*c:
			return 0.5
		else:
			return (a*d - math.sqrt(a*b*c*d))/((a*d)-(b*c))

def survivorFixed(k,M,n,N):
    """
    Evaluates the Hypergeometric log-proability mass function for drawing k certain residues from a sample of N
    given M objects of which n are of this same type

    For each Tuminello possible Node pair, this is calculated for:
    
    tummData[0]-1,len(msa),tummData[1],tummData[2]

    = (wProj - 1, full ali seq number,Dj, Di)

    = (Res 2 given Res 1 seq number - 1, full ali seq number,Res 2 seq number, Res 1 seq number)
    =  k,                              M,                  n,               N

    So this function:

    Evaluates the Hypergeometric log-proability mass function for drawing Res 2 given Res 1 seq number from a sample of Res 1 seq number
    given the full ali seq number objects of which Res 2 seq number are of Res 2 type

    Or observed from a sample given observed in full ali

    """
	quant, tot, good, draw = k, M, n, N
	k2 = np.arange(quant + 1, draw + 1)
	logpmf = hypergeom._logpmf(k2, tot, good, draw)   # Evaluate the log-pmf instead
	logsf = logsumexp(logpmf)
	return logsf


def survivorFixedWithAntiCorr(k,M,n,N):
    """
    Evaluates the Hypergeometric log-proability mass function for drawing k certain residues from a sample of N
    given M objects of which n are of this same type

    For each Tuminello possible Node pair, this is calculated for:
    
    tummData[0]-1, len(msa), tummData[1], tummData[2]

    = (wProj - 1, full ali seq number,Dj, Di)

    = (Res 2 given Res 1 seq number - 1, full ali seq number,Res 2 seq number, Res 1 seq number)
    =  k,                                M,                  n,               N

    So this function:

    Evaluates the Hypergeometric log-proability mass function for drawing Res 2 given Res 1 seq number from a sample of Res 1 seq number
    given the full ali seq number objects of which Res 2 seq number are of Res 2 type

    Or observed from a sample given observed in full ali
    #NEW: With AntiCorre
    """
	quant, tot, good, draw = k, M, n, N
    k2 = None
    if k > (N * (n/M)):
        k2 = np.arange(quant + 1, draw + 1)
    else:
        k2 = np.arange(0, quant + 1)
	# k2 = np.arange(quant + 1, draw + 1)
	logpmf = hypergeom._logpmf(k2, tot, good, draw)   # Evaluate the log-pmf instead
	logsf = logsumexp(logpmf)
	return logsf

def getTumminelloTuple(res1,res2):
    """
    Calculate Tumminello edge correlation for two residue Node objects (res1, res2)
    """
	aa1_list = conan_constants.setsAA[res1.id[0]]
	i1 = res1.position[0]-1
	aa2_list = conan_constants.setsAA[res2.id[0]]
	i2 = res2.position[0]-1
	wProj = 0
	Di = 0
	Dj = 0
	a = 0.0
	b = 0.0
	c = 0.0
	d = 0.0

    #Counts are now done in the whole MSA for:
    # res1 (Di), res2 (Dj)
    # res2 if res1 present (wProj and a)
    # res1 if res2 not present (b)
    # res2 if res1 not present (c)
    # neither res1 nor res2 present (d)
	for sequence in msa.values():
		if sequence[i1].upper() in aa1_list:
			Di +=1.0
			if sequence[i2].upper() in aa2_list:
				Dj += 1
				wProj += 1
				a += 1.0
			else:
				b += 1.0
		elif sequence[i2].upper() in aa2_list:
			Dj += 1
			c += 1.0
		else:
			d += 1.0

    #jaccard's coefficient calculated
	jc = a/(a+b+c)

	if Di > Dj:
        #for possible correlation, returns frequencies and jaccard
		return ((wProj,Di,Dj),jc)
	else:
        #for possible anticorrelations or no correlations, also returns frequencies and jaccard
		return ((wProj,Dj,Di),jc)