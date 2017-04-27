from __future__ import division
import math
##############################################################################
# SFS FROM FREQUENCIES
# This function takes a vector of frequencies (a python list) and turns them into a SFS. 
# You must provide the length of the region from which the sites were taken and the number of alleles sampled
# Nice and general
##############################################################################

def SFS_from_frequencies(frequencies, length,N):
	SFS = [0]*(N+1)
	for i in frequencies:
		if i > N:
			print "SFS_from_frequencies: Error in your frequencies vector: One of the values is greater than the number of individuals\nThe offending value is: " + str(i) +" and the sample is "+str(N)
			return
		SFS[i] += 1
	SFS[0] = length - len(frequencies)
	if sum(SFS) < length:
		print"SFS_from_frequencies: Error in your frequencies vector: Fewer items in the SFS than the length of the region"
		return
	if sum(SFS) > length:
		print"SFS_from_frequencies: Error in your frequencies vector: More items in the SFS than the length of the region"
		return
	return SFS


def SFS_from_all_frequencies(frequencies,N):
	SFS = [0]*(N+1)
	length = len(frequencies)
	for i in frequencies:
		if i > N:
			print "SFS_from_frequencies: Error in your frequencies vector: One of the values is greater than the number of individuals\nThe offending value is: " + str(i) +" and the sample is "+str(N)
			return
		if i > 0 and i <= N:
			SFS[i] += 1
			
	SFS[0] = length - sum(SFS)
	if sum(SFS) < length:
		print"SFS_from_frequencies: Error in your frequencies vector: Fewer items in the SFS than the length of the region"
		return
	if sum(SFS) > length:
		print"SFS_from_frequencies: Error in your frequencies vector: More items in the SFS than the length of the region"
		return
	return SFS

##############################################################################
# MERGE SFS
# This function takes a two SFS and sums their contents. 
# They must be of equal size (i.e. same N)
# Use this to get SFS for corresponding bins within a slim run
##############################################################################
def merge_SFS(a,b):
	import operator
	if len(a) != len(b):
		print "merge_SFS: The length of the two SFS are different"
		return		
	c = map(operator.add,a,b)
	return c


##############################################################################
##############################################################################
##############################################################################
#
#		The following are a bunch of statistics that can be extracted 
#		from the SFS
#
##############################################################################
##############################################################################
##############################################################################

# Pi
# This function returns the estimate of theta_pi from a SFS 
# 
##############################################################################


def pi(SFS,per_site = True):
	if sum(SFS) ==0: 
		return -99
	N = len(SFS)-1
	binom = (N * (N -1))/2
	pi = sum([(1.0*i*(N-i)*(SFS[i]))/(binom) for i in xrange(N) if i != 0])
	if per_site == True:
		return pi/sum(SFS)
	else:
		return pi

##############################################################################
# Xsi
# This function returns Peter's Xsi statistic, a stataistic that could 
# potentially be diagnostic for adaptive evolution 
# 
##############################################################################

def xsi(SFS):
	if len(SFS) ==0: return -99
#        print SFS
	N = len(SFS)-1
	total = sum(SFS[1:-1])
	xsi_numerator = sum([float(i*i*SFS[i])/total for i in xrange(N) if i > N/2 and i !=N])
	xsi_denominator = sum([float(i*SFS[i])/total for i in xrange(N) if i > N/2 and i !=N])
	if xsi_denominator == 0:
		return -99
	else:
		return (xsi_numerator/xsi_denominator)/N

##############################################################################
# Pi2
# This function returns the estimate of theta_pi from a SFS 
# 
##############################################################################


def pi2(SFS):
	if len(SFS) ==0: return -99
	N = len(SFS)-1
	binom = (N * (N -1))/2
	pi2 = sum([(1.0*i*(N-i)*(SFS[i]))/binom for i in xrange(N) if i >= N/2 and  i != N ] )
	return (pi2/sum(SFS))/pi(SFS)

##############################################################################
# fwh
# Returns Fay and Wu's estiamte of the population mutation parameter [theta] 
# 
##############################################################################


def fwh(SFS,per_site = True):
	if len(SFS) ==0: return -99
	N = len(SFS)-1
	binom = (N * (N -1))/2
	fwh = sum([(1.0*i*i*(SFS[i]))/binom for i in xrange(N) if i != 0 or  i != N] )
	if per_site == True:
		return fwh/sum(SFS)
	else:
		return fwh

##############################################################################
# Theta_W
# Returns Watterson's estimator of the mutation rate
# 
##############################################################################

def theta_W(SFS,per_site=True):
	if len(SFS) ==0: return -99
	N = len(SFS)-1
	S = sum(SFS[1:-1]) ## this slice takes the interior of the SFS, gets S
	harmonic = sum(1.0/d for d in xrange(1, N-1,1))
	if per_site == True:
		return float(S)/(harmonic*sum(SFS))
	else:
		return float(S)/(harmonic)

##############################################################################
# tajima
# Returns Tajima's D (198?) The estimate of the skew of the SFS
# 	~~Based on Rob's script, mine had a cryptic typo
##############################################################################

def tajima(SFS):
	#if len(SFS) <2: return None
	th_pi = pi(SFS,per_site=False)
	N = len(SFS)-1
	S = sum(SFS[1:-1])
	length = sum(SFS)
	a1 = sum([1.0/i for i in range(1, N)])
	a2 = sum([1.0/(i**2) for i in range(1, N)])
	b1 = float(N+1)/(3*(N-1))
	b2 = float(2 * ( (N**2) + N + 3 )) / (9*N*(N-1))## Tajima 1989 Equation 9
	c1 = b1 - 1.0/a1
	c2 = b2 - float(N+2)/(a1 * N) + float(a2)/(a1 ** 2)
	e1 = float(c1) / a1
	e2 = float(c2) / ( (a1**2) + a2 )
	if ((e1 * S )+ ((e2 * S) * (S - 1))) == 0.0:   # Tajima 1989 Equation 38
		return -99
	else: 
		D = (float(th_pi - (float(S)/(a1)))/math.sqrt((e1 * S )+ ((e2 * S) * (S - 1) )))
		return D


