import shared_globals
from decimal import Decimal
import numpy as np
import math
import conan_constants

def lnfact(x):
	"""
	Returns log in base e of Factorial
	"""
	result = np.int(10)
	if x==0 or x==1:
		return 0.0
	else:
		for i in range(1,x+1):
			result = result*i
	return math.log(result)

def stirling(x):
	"""
	Used as base for Stirling's Normal Approximation for binomial distribution

	For x below 9 returns log in base e of Factorial
	For x equal or above 9 returns complicated number

	"""
	if x<=8:
		return lnfact(x)
	else:
		return ((x+0.5)*math.log(x))-x+0.918938533205

def lnbdf(N,nx,px):
	"""
	Returns p-value for binomial distribution using Stirling's Normal Approximation for binomial distribution
	"""
	if N==nx:
		return N*math.log(px)
	if nx==0:
		return N*math.log(1.0-px)
	else:
		return stirling(N)-stirling(nx)-stirling(N-nx)+nx*math.log(px)+(N-nx)*math.log(1-px)

def eto10(x):
	"""
	Converts base e number to base 10 (decimal)
	"""
	return x/2.30258509299

def cbd_tietjen(N,n,freq,right):
	"""
	Cumulative Binomial Distribution calculus using Tietjen algorithm
	https://www.tandfonline.com/doi/abs/10.1080/00031305.1994.10476043
	Results uses the Decimal class and thus small numbers are only lower-limit bounded by memory
	"""
	a = math.pow(1-freq,N)
	if not right and n== 0:
		return a
	if right and n == 0:
		return 1
	suma = Decimal(0.0)
	if not right:
		suma = Decimal(a)
	u = freq/(1-freq)
	for i in range(2,int(N+2)):
		a = float(a*u*(N+2-i))/float(i-1)
		if (not right and i <= n+1) or (right and i-1 >= n):
			suma += Decimal(a)
	return suma

def cbd(N,n,freq,right):
	"""
	Calculates Cumulative Binomial Distribution for each individual case by different algorithms (see comments for details)

	Calculus can be done either using Tietjen algorithm or Stirling's Normal Approximation
	"""
	#integer conversion of counts
	n = int(n)
	N = int(N)
	#if 1-freq to the power of sample space N is not 0
	# (smaller frequency in relatively low N counts or higher frequency in low/high N)
	if math.pow(1-freq,N) != 0:
		#calculate value with cumulative binomial distribution tietjen
		val = cbd_tietjen(N,n,freq,right)
		if val == 0.0:
			#if value is really small, return this very small number
			return 5e-324
		else:
			#if value is not terribly small, convert result to log10 base p-value score
			return math.ceil(math.log10(cbd_tietjen(N,n,freq,right)))
	#if we are working with either really high N counts or quite low residue frequencies
	suma = 0
	result = 0
	ps = []
	#P-value is calculated by Stirling's Normal Approximation for binomial distribution in succession for each count
	minP = eto10(lnbdf(N,n,freq))
	#right means AUC from n until N
	if right:
		for i in range(n,N+1):
			val = eto10(lnbdf(N,i,freq))
			ps.append(val)
			if val > minP:
				minP = val
		for val in ps:
			suma += math.pow(10,val)-math.floor(minP)
	#else means AUC from 0 until n
	else:
		for i in range(0,n+1):
			val = eto10(lnbdf(N,i,freq))
			ps.append(val)
			if val > minP:
				minP = val
		for val in ps:
			suma += math.pow(10,val)-math.floor(minP)
	return math.floor(math.log(suma))+math.floor(minP)

def DRCN(res1,res2):
	"""
	DRCN Correlation calculation method for two given Node objects (res1, res2)
	"""
	# global msa

	aa1_list = conan_constants.setsAA[res1.id[0]]
	i1 = res1.position[0]-1
	aa2_list = conan_constants.setsAA[res2.id[0]]
	i2 = res2.position[0]-1
	Na = 0.0
	Nb = 0.0
	Nba = 0.0
	a = 0.0
	b = 0.0
	c = 0.0
	d = 0.0
	n = float(len(shared_globals.msa))

	for sequence in shared_globals.msa.values():
		if sequence[i1].upper() in aa1_list:
			Na +=1.0
			if sequence[i2].upper() in aa2_list:
				Nb += 1.0
				Nba += 1.0
				a += 1.0
			else:
				b += 1.0
		elif sequence[i2].upper() in aa2_list:
			Nb += 1.0
			c += 1.0
		else:
			d += 1.0

	jc = a/(a+b+c)

	if Nba > (Na*(Nb/n)):
		freq1 = float(Nb)/float(n)
		freq2 = float(Na)/float(n)
		pv1 = cbd(Na,Nba,freq1,True)*-1
		pv2 = cbd(Nb,Nba,freq2,True)*-1

		if pv1 > pv2:
			return (pv1,jc)
		else:
			return (pv2,jc)
	else:
		#Anti-Correlation
		return -1.0

def DRCNWithAntiCorr(res1,res2):
	"""
	DRCN Correlation calculation method for two given Node objects (res1, res2)
	NEW: AntiCorrelations added
	"""
	# global msa

	#conversion from residue id to AA Type or Column Property
	aa1_list = conan_constants.setsAA[res1.id[0]]
	i1 = res1.position[0]-1
	aa2_list = conan_constants.setsAA[res2.id[0]]
	i2 = res2.position[0]-1
	
	#Count initialization
	Na = 0.0
	Nb = 0.0
	Nba = 0.0
	a = 0.0
	b = 0.0
	c = 0.0
	d = 0.0
	n = float(len(shared_globals.msa))

	#Counts are now done in the whole MSA for:
	# res1 (Na), res2 (Nb)
	# res2 if res1 present (Nba and a)
	# res1 if res2 not present (b)
	# res2 if res1 not present (c)
	# neither res1 nor res2 present (d)
	for sequence in shared_globals.msa.values():
		if sequence[i1].upper() in aa1_list:
			Na +=1.0
			if sequence[i2].upper() in aa2_list:
				Nb += 1.0
				Nba += 1.0
				a += 1.0
			else:
				b += 1.0
		elif sequence[i2].upper() in aa2_list:
			Nb += 1.0
			c += 1.0
		else:
			d += 1.0

	#jaccard's coefficient calculated
	jc = a/(a+b+c)

	#if res2 given res1 count is BIGGER than the Expected res1 count multiplied by res2 frequency
	if Nba > (Na*(Nb/n)):
		#freqs are calculated
		freq1 = float(Nb)/float(n)
		freq2 = float(Na)/float(n)
		# print(Na,Nba, n, freq1, freq2)

		#cumulative binomials are calculated for right
		pv1 = cbd(Na,Nba,freq1,True)*-1
		pv2 = cbd(Nb,Nba,freq2,True)*-1
		# print("here 4")

		if pv1 > pv2:
			return (pv1,jc, freq1-freq2)
		else:
			return (pv2,jc, freq2-freq1)
	#if res2 given res1 count is smaller than the Expected res1 count multiplied by res2 frequency
	# else:
		# return -1.0
		# return -1.0, jc
	elif Nba < (Na*(Nb/n)):
		#freqs are calculated
		freq1 = float(Nb)/float(n)
		freq2 = float(Na)/float(n)
		# print("anti")
		# print(Na,Nba, n, freq1, freq2)

		#cumulative binomials are calculated for left
		pv1 = cbd(Na,Nba,freq1,False)#*-1
		pv2 = cbd(Nb,Nba,freq2,False)#*-1
		if pv1 < pv2:
			return (pv1,jc, freq1-freq2)
		else:
			return (pv2,jc, freq2-freq1)
	#if res2 given res1 count is equal to the Expected res1 count multiplied by res2 frequency
	else:
		return 0.0, jc, 0.0