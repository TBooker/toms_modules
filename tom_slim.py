##############################################################################
# SLIM_READER
# For large files containing many slims which can be very large (~1Mb per sim)
# It is pretty inefficient to read these in as a whole file then split on 
# some value. It is better to iterate through the file one line at a time and
# choose certain lines to give out (that's what the yield call does)
##############################################################################

## write a text parser that allows for effient reading of a large slim output collection
def slim_reader(input_file):
	with open(input_file,"r") as FILE:
		slim_output = ''
		terminator = 'None\n'
		for line in FILE:
			if line.startswith('#INPUT'):
				slim_output += line
			elif not line.startswith(terminator):
				slim_output += line
			elif line.startswith(terminator):
				yield slim_output
				slim_output = ' '
			
		#yield slim_output		
import gzip
def slim_reader_gzip(input_file):
	with gzip.open(input_file,'r') as FILE:
		slim_output = ''
		terminator = 'None\n'
		for line in FILE:
			if line.startswith('#INPUT'):
				slim_output += line
			elif not line.startswith(terminator):
				slim_output += line
			elif line.startswith(terminator):
				yield slim_output
				slim_output = ' '
			
		#yield slim_output	


##############################################################################
##############################################################################
#
# SLIM > class
# I decided to make a slim file into a class oject so that I would only have 
# to calculate and extract certain things once. It is maybe a little bit more
# memory intensive, but th processing savings make up for this I think.  
# ONLY HANDLES A SINGLE TIME POINT
##############################################################################
##############################################################################
######################################################################
###
###	This class should provide a (relatively) flexible way to analyse slim files
###	By supplying a slim file in either string or list format to the class
### 	you get an object that has all of the variables and items of interest from the 
###	slim file. By doing it this way, you only need to iterate through the slim file
###	once and 
###
###		The slim class has several methods:
###
###	slim.name 		name of the file this comes from << TO DO
###	slim.N			Number of individuals in the whole population
###	slim.length		The length of the simulated chromosome
###	slim.sampleN		Number of indiviuduals in the sample (if not random sample
###					this i qual to slim.N
###	slim.theta		theoretical population mutation rate (4Neu)
###	slim.sites_dict		a dictionary of selected/neutral sites for each site type
###					the keys are "selected" and "neutral"
###					these access lists of the selected sites
###					for each site type
###	slim.regions		a list of region types
###	slim.organs		chomosome organistation
###	slim.DFE 		a list of the differnt DFE usedin the SLiMulation
###	slim.mut_rate		the mutations/site/generation rate	
###	slim.recom_rate		the recombination rate	as a list of recombination
###					intervals. Could add in a function to get a 
###					weighted mean of the recombination rate for 
###					region
###	slim.rho		population recombination rate (4Ner)
###	self.generations	number of generations
###	self.mutations		the mutations as a list of lists
###	self.fixed		The fixed mutations as a list of lists
###	
### ** could add in demography later if Sims become more complicated
###
######
######
######		The class has a couple of convienience functions...
###

###	slim.sanity_test(threshold=0.0001)This function returns the string "good" if all 
###					of the mutations in the simulation are behaving
###					 themselves. If not, it returns a None 		
###	slim.sites_dict()		This function returns a dictionary with keys: "selected"
###					and "neutral". These are lists with all the sites 
###					of either type in them
###
###	slim.organ_lengths()		Gives a dictionary of the combined lengths of the genomic
###					elements in a particular 
###
###	slim.organs_mutations()		Gives a dictionary of the mutations in a particular
###					genomic element, uses the SLiM codes (e.g. g0, g1...)
###
######################################################################

from tom import in_range

"""This is a class that allows you to access all the elements of a slim run easily
For example, if you want to look at how chromosomal structure and genetic diversity interact
you can get the information easily. 
There are more options on the way...
Fixed option allows you to specify whether ther are fixed mutations in your output"""
class slim:
	def __init__(self,slimput,fixed=False,give_genomes = False): 
		if type(slimput) == str:
			slimmers = slimput.split("\n")	
		elif type(slimput) == list:
			slimmers = slimput
		if len(slimmers) < 10:
			print("Error 1 for class slim: SLiM input is too short")
			self.sanity = None
			return None

		self.sanity = "sane"
		self.all = slimmers
		self.name = slimmers[slimmers.index("#INPUT PARAMETER FILE")+1].strip("\n")
		self.N = int(slimmers[slimmers.index('#DEMOGRAPHY AND STRUCTURE')+1].split(" ")[3])  ## return the number of individuals in the shol_simulation
		output_lines = [i for i in slimmers if i.startswith('#OUT:') or i.startswith("G")]
		if len(output_lines) ==0:
			self.output = False
			print "File not output"
			return None
		else:
			self.output = True
		random_line = [i for  i in output_lines if "R" in i][0]
		self.sampleN = int(random_line.split(" ")[4])
		fixed_line = slimmers.index([i for  i in output_lines if "F" in i][0])
		genome_line = slimmers.index([i for i in output_lines if "G" in i][0])
		if fixed:		
			if len(fixed_line) == 0:
				fixed_muts = False
			else:
				fixed_muts = True

		self.mut_rate = float(slimmers[slimmers.index('#MUTATION RATE')+1])
		self.theta = float(self.N) * self.mut_rate * 4.0
		self.generations = int(slimmers[slimmers.index('#GENERATIONS')+1])
	#######################################################################
		recomb =[]
		recom_index = slimmers.index("#RECOMBINATION RATE")
		for i in xrange(len(slimmers)):
			ind = i + 1
			if slimmers[recom_index+ind].startswith("#"):
				break
			else:
				recomb.append(slimmers[recom_index+ind].strip("\n").split(" ")) 
		self.recomb_intervals = recomb
	#######################################################################
		DFE_index =  slimmers.index('#MUTATION TYPES')  
### The DFE will be returned as a list as there may be more than one
		the_DFEs = []
		for line in slimmers[DFE_index+1:]:
			if line.startswith("#"):break
			else: the_DFEs.append(line.strip("\n").split(" "))
		self.DFE= the_DFEs ## Now a list of the DFEs

	#######################################################################
		organ_index =  slimmers.index('#CHROMOSOME ORGANIZATION')  	
		
		the_organs = []
		for line in slimmers[organ_index+1:]:
			if line.startswith("#"):break
			else: the_organs.append(line.strip("\n").split(" "))
		self.organs= the_organs ## Now a list of the regions

	#######################################################################			
		region_index =  slimmers.index('#GENOMIC ELEMENT TYPES')  	
		the_regions = []
		for line in slimmers[region_index+1:]:
			if line.startswith("#"):break
			else: the_regions.append(line.strip("\n").split(" "))
		self.regions= the_regions ## Now a list of the region
	#######################################################################
		self.length = sum([(int(l[2])-int(l[1]))+1 for l in self.organs])
	#######################################################################	
		if random_line:
			mutations_index = slimmers.index(random_line)+1
			mutations_raw = [i.split(" ") for i in slimmers if slimmers.index(i) > mutations_index and slimmers.index(i) < genome_line and len(i.split(" ")) == 8] ## this relies on the structure of the output HEAVILY!
## Thats a doozy, but it gets all the mutations you need from the populations
			self.mutations = mutations_raw
		else:
			self.mutations = None
		
		if fixed:
			if fixed_muts:
				self.fixed = [i.split(" ") for i in slimmers[:] if len(i.split(" ")) == 8]
			elif not fixed_muts:
				self.fixed = []
		if give_genomes:
			self.genomes = [i for i in slimmers[genome_line+1:fixed_line]]
		#	mutations_raw = [i.split(" ") for i in slimmers if slimmers.index(i) > mutations_index and slimmers.index(i) < slimmers.index() and len(i.split(" ")) == 8]
		else:
			self.genomes = None
################################################################################		
	#################################################
### This function tests whether the output is sensible or not
	def slim_sanity_test(self, threshold = "0.0001"):
		if self.sanity == None:
			print "Error 1 for class slim: SLiM input is too short"
			return None
		else:
			mut_info = [[float(line[7])/float(self.N),line[3],line[6]] for line in self.mutations if float(line[3]) != 0]  # get the allele frequency of all mutations
			if len(mut_info) == 0:
				print "No non-neutral mutations in the population"
				return "good"
			for i in mut_info: 	#|||| i is frequency |||| j is selection ||||
				q = i[0]
	
				s = i[1]
				gen = i[2] 

				if float(s) <= -2 and float(q) > threshold and int(gen) < int(self.generations):
					#print(i)
					print "Error 3: SLiM run seems to have instances of the weird selection/frequecny discrepancy"			
					return None
					
				else:
					return "good"
################################################################################
## SITES DICT PSEUDO CODE...
#
## This function is to take a slim and return a dictionary of sites.
# the keys of which are selected and neutral. The neutral part has the sites that are neutral in the whole simulated chromosoem
# for element type:
#	if element type is not the interelement neutral site type:
#		open list
#		for element on the chromosome
#			if this element is a selected type
#				add the sites in this region to the list
#		the list is now populated with all the sites of a partivular element type
#		add this to the selected type list with a reference to the particular elemetn type
# make the selected element of the dictionary correspond to the list of elements and their positions
#
#	repeat for neutral sites
#
#
#

## Not a very good way of doing this, but good for the HRI scripts...
	def sites_dict(self):
		sites_dict = {}
		selected_sites =[]
		for j in self.regions:
			if j[0] == "g0":continue
			selected_sites_x = []

			for i in self.organs:
				if i[0] == j[0]:
					selected_sites_x += range(int(i[1]),int(i[2])+1)
			selected_sites.append([j,selected_sites_x])
		sites_dict["selected"] = selected_sites		
		
		neutral_sites = []		
		for i in self.organs:
			if i[0] == "g0":  # BREAKABLE, but efficient, ALWAYS have neutral sites as g0
				neutral_sites += range(int(i[1]),int(i[2])+1)
		sites_dict["neutral"] = ["neutral",neutral_sites]
		return sites_dict
################################################################################
	def organ_lengths(self):
		organ_lengths = {}
		for org in self.organs:
			if org[0] not in organ_lengths.keys():
				organ_lengths[org[0]]=int(org[2])-int(org[1])+1
			else:
				organ_lengths[org[0]]+=int(org[2])-int(org[1])+1
		return organ_lengths
################################################################################
	
	def organ_mutations(self):
		org_dict = {}
		for mut in self.mutations:
		#	print mut
			for org in self.organs:
				if in_range(int(mut[2]),[int(org[1]),int(org[2])]):
					if org[0] not in org_dict.keys():
						org_dict[org[0]] = [mut]
					else:
						org_dict[org[0]].append(mut)
		return org_dict
################################################################################
	def mutations_dict(self):
		mut_dict = {}
		for mut in self.mutations:
			mut_dict[int(mut[0])] = int(mut[2]) 
			
		return mut_dict

	def genome_dict(self):
		genome_dict={}
		for i in self.genomes:
			x = i.split(" ")
			genome_dict[x[0]] = x[1:]
		return genome_dict
