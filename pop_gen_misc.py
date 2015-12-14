##############################################################################
# R_ZETA
# This function returns Riemann's Zeta function for a given shape parameter 
# Could be sped up if it were implemented in Cython, but since it is usede so 
# infrequently I don't see the point
##############################################################################

def r_zeta(s,max_iter=1000000,tolerance = 0.0000000001):
	zeta_list = []
	for j in xrange(max_iter):
		if j == 0:continue
		if len(zeta_list) ==0:		
			zeta = j**-(s)
			zeta_list.append(zeta)
		else:
			zeta = j**-(s)
			diff = zeta_list[-1] - zeta
			if diff < tolerance:
				break
			else:
				zeta_list.append(zeta)
	print "converged on a value of "+str(sum(zeta_list))+" after "+str(j) + " iterations"
	return sum(zeta_list)			 


##############################################################################
# jukes_cantor
# This function applies the jukes-cantor correction for multiple hits 
#
# COULDN'T THINK OF A BETTER PLACE FOR THIS
##############################################################################
from math import log
def jukes_cantor(raw):
	return -0.75*(log(1-(4*float(raw)/3)))




##############################################################################
# WELCH_TRANSFORM
# This function applies the transformation of a gamma DFE inferred by DFE-alpha 
# into the one recommended by BC
# 
##############################################################################

def welch_transform(sh,o_d,zeta):
	return pow(o_d/(sh**(sh+1)*zeta),-1.0*sh)
