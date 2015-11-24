#### Make a module of functions that I can call at other points in time
#### Basically, generic things that I use a lot in the slim analysis
#### it is important that all of the fnctions in here are as flexible and efficient as possible

import os,glob,sys, math
import numpy as np

##############################################################################
# MERGE FILES
# This function takes a file identifier and finds all files that match it. 
# Then it combines them all into a single file 
# Good for parallelised results and for making memory savings!
# UPDATE:
#	 Allows you to specify whether to remove redundant headers
#	 Adds a column telling you which file the rows came from
##############################################################################

def merge_files(ident,merge_name,clean="F",header = "F"):
	files = glob.glob( ident+"*" )
	with open( merge_name, 'w' ) as result:
	    if header =="T":
		tick ="Y"
	    for file_ in files:
		if header.upper() == "T":
			line_num =1
		elif header.upper() != "T":
			line_num =2
	        for line in open( file_, 'r' ):
		    if line_num == 1:
			if tick !="N":
				result.write("file,"+line)
				tick ="N"
			line_num+=1
			continue
		    else:
	        	result.write(file_ + "," +  line )
			line_num +=1			
	if clean == "F": return
	else:	
		os.system("rm " + ident +"*")
		return
##############################################################################
# BRACE
# This function just asks you to hit enter or to quit the program. Good for 
# developing scripts. Tres useful
##############################################################################
def brace():
	x = raw_input("Hit enter to pass or enter any key to quit: ")
	if not x:
		print("PASS")
		pass
	elif x:
		print("EXIT")
		sys.exit()

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
# OVERLAP
# this function takes two sets of [start, stop] and tells you whether they
# overlap. If they do, it returns the limits of the overlap.
##############################################################################


def overlap(a,b):
## If non_overlap == True, returns the non_overlapping section
## Otherwise it returns the overlapped sections
## If non_overlap == False and the two sections are identical, returns None
	if a[1] >= b[1] and a[0] <= b[1]:
		x = [min([a[0],b[0]]),max([a[1],b[1]])]
	elif b[1] >= a[0] and b[0] <= a[1]:
		x = [min([a[0],b[0]]),max([a[1],b[1]])]
	elif a[0] <= b[0] and b[1] <= a[1] or b[0] <= a[0] and a[1] <= b[1]:
		x = [min([a[0],b[0]]),max([a[1],b[1]])]	
	elif a[0] == b[0]: 
		x = a
	else:
		x = None
	return x

##############################################################################
# BOOTSTRAP SAMPLE
# If you give this function a number, corresponding to the entries you want to bootstrap
# It returns a list of samples. it Uses numpy
##############################################################################


def bootstrap_sample(number_of_entries):
	sample_array = np.random.rand(number_of_entries)*number_of_entries
	sample_list = map(int,sample_array)
	return sample_list



##############################################################################
# MEAN_SE
# I wrote this for getting the mean of stats I extracted from Slim output
# but it should be flexible enogh to calculate mean and standard error for
# any list of numbers. The index option allows you to pass a list of lists and 
# to choose an index of the internal lists
##############################################################################


def mean_se(stats_raw,index = "NA"):
	if index == "NA":
		stats = stats_raw 
	elif type(index) == int:
		stats = [stat[index] for stat in stats_raw if stat[index] != "NA" and stat[index] !=None]
	else:
		print("The index was supplied innapropriately")
		return
	if len(stats) > 0 :
		n = len(stats)
		mean = sum([float(stat) for stat in stats if stat])/n
		sum_squares= sum([(float(stat)-mean)**2 for stat in stats])
		sd = math.sqrt(sum_squares/(n-1))  # sampling variance should be (n-1)
		se= sd/math.sqrt(n)
		return mean,se,n
	else:
		mean="NA"
		se="NA"
		n ="NA"
		return mean,sw,n

##############################################################################
# ALL_SAME
# This little function just returns TRUE if all the items in a list are the 
# same. It can also be acheived by converting the list to a set and then testing 
# the length 
##############################################################################

def all_same(items):
    return all(x == items[0] for x in items)

##############################################################################
# IN_RANGE
# This function returns True if a given int falls within a given range 
##############################################################################

def in_range(point,range):
	if point >= range[0] and point <= range[1]:
		return True
	else:
		return False
