###
### My module of useful VCF functions. These are, however, not as robust as those found in PyVCF. 
### They are however, faster
###



##############################################################################
# VCF_LINE
# This is a class that converts a VCF line as a string or as a list and get
# information.
#
# The object is similar to a PyVCF line object, but is a bit more fragile than
# PyVCF as it does not look at the header, rather, it assumes that everything 
# is in a standard order.
##############################################################################
class VCF_line:
	""" This is my VCF_line class, it is basically just a mini version of the pyVCF record object"""
	""" It allows you to access a VCF as a Tabix object instead of having to make a big file object"""  

	def __init__(self,line):
		if type(line) == str:
			w = line.split("\t")
		elif type(line) == list:
			w = line
		else:
			print "unsuitable line given to the parse_VCF_line function"
			return None
		self.pos = w[1]
		self.chrom=w[0]
		self.qual = w[5]
		
		info_dict ={}
		filter_dict = {}
		for i in w[7].split(";"):
			filt = i.split("=")	
			if len(filt) == 1:
				filter_dict[filt[0]]= None
			else:
				filter_dict[filt[0]]=filt[1]
		self.filters = filter_dict
		if "INDEL" in filter_dict.keys():
			self.indel = True
		else:
			self.indel = False
		allele_dict={}
		allele_dict["0"] = w[3]
		info_fields = w[8].split(":")
		if w[4] !=".":
			self.variant = True
			for q,l in zip(w[4].split(","),xrange(1,len(w[4].split(",")))):
				if q == "X":continue
				allele_dict[str(l)] = q
		else: self.variant=False
		self.alleles = allele_dict
		for o in info_fields:
			info_dict[o] =[]
		self.nsamples = len(w[9:])
		for r in w[9:]:
			info = r.split(":")
			for g in info_dict.keys():
				info_dict[g].append(info[info_fields.index(g)])
		self.info = info_dict
	def total_depth(self):
		return sum([float(dp) for dp in vcf_line.info["DP"]])

##############################################################################
# VCF_allele_frequencies
# Using a VCF_line object as defined above, the function returns a list of 
# alleles (e.g. [A,G,G,G,A,A,A,G]
# As it is now, you can specify GQ and DP to filter sites on 
#
##############################################################################
def VCF_allele_frequencies(vcf_line,ploid = 2,minGQ=0,minDP=0):
	""" This function will return a dict of the alleles and their respective frequencies"""
	allele_list = []
	if not vcf_line.variant:
		for dp in vcf_line.info["DP"]:
			if float(dp) < minDP: continue
			for i in xrange(ploid):
				allele_list.append(vcf_line.alleles["0"].upper())
	else:
		for geno,gq,dp in zip(vcf_line.info["GT"],vcf_line.info["GQ"],vcf_line.info["DP"]):
			if float(gq) < minGQ or float(dp) < minDP :continue
			allele_list.append(vcf_line.alleles[geno.split("/")[0]].upper())		
			allele_list.append(vcf_line.alleles[geno.split("/")[1]].upper())		
	ref_allele = vcf_line.alleles["0"]
	if ref_allele not in allele_list: 
		return 
	return allele_list

##############################################################################
# alleles_2_PKformat
# This converts a list of allele frequencies into a string that is accepted by 
# PK's sfs inference programs
# 
##############################################################################
def alleles_2_PKformat(allele_list):
	A = str(allele_list.count("A"))
	C = str(allele_list.count("C"))
	G = str(allele_list.count("G"))
	T = str(allele_list.count("T"))
	return " ".join([A,C,G,T])


##############################################################################
# CONVERT_ONE_LINE
#
# Takes a tab-separated string of allele frequencies (or list of allele frequencies)
# and converts in to the PK format
#
# Pretty specialised and not really transferable to other projects
##############################################################################
def convert_one_line(line):
	if type(line) == str:
		line.split("\t")
	elif type(line) == list:
		pass		
	else:
		print "incorrect object type provided"
		return
	out_string = ""
	for s in line:
		if s == ".,.,.,.":
			return
		x = s.split(",")
		if len(set(x)) >3:  # Remove any site with > 2 alleles at a given position
			return
		out_string += " ".join(x)+" "
	return out_string

##############################################################################
# DIVERGENT
#
# Takes a list of PK format allele frequencies and returns 1 if the outgroup is
# divergent, 0 if not. None for a triallelic site 
#
##############################################################################
def divergent(ingroup,outgroup,out_alleles = 1):
	if len(set(ingroup)) >3:  # Remove any site with > 2 alleles at a given position
		return
	if 20 in ingroup:
		if ingroup.index(20) != outgroup.index(out_alleles):
			return 1
		else:
			return 0
	if 20 not in ingroup: ## polymorphic
### Takes a frequency weighted allele
		items = [i for i in ingroup if i !=0]
		choosing_list = [items[0]]*items[0] + [items[0]]*items[0] 
		allele_chosen = random.choice(choosing_list)
		if ingroup.index(allele_chosen) != outgroup.index(out_alleles):
			return 1
		else:
			return 0



