import argparse, sys, pysam, gzip, random, pickle
from collections import Counter
import site_frequency_spectrum as SFS_tools

def convert_freq_alleles(alleles):
# Returns a list of alleles and their frequencies
	allele_dict = {0:'A',1:'C',2:'G',3:'T'}
	z = map(int,alleles.split(','))
	conv = [[allele_dict[i],z[i]] for i in range(4) if z[i] != 0]
	return conv

def random_allele(alleles):
	listy = []
	for i in alleles:
		listy += i[0]*i[1]
	return random.choice(listy)

def polymorphic(alleles):
	if len(alleles) > 1:
		return True
	else:
		return False

class freq_line:
	def __init__(self,split_line): 
		self.raw = '\t'.join(split_line)
		self.chrom = split_line[0]
		self.pos = int(split_line[1])
		self.depth = int(split_line[2])
		self.HWE = split_line[3]
		self.cast_cpg = split_line[4]
		self.fam_cpg = split_line[6]
		self.rat_cpg = split_line[8]

		if self.cast_cpg == '.':
			self.cast_called = False
			self.cast_alleles = None
		else:
			self.cast_called = True
			self.cast_alleles = convert_freq_alleles(split_line[5])
			
		if self.fam_cpg == '.':
			self.fam_called = False
			self.fam_alleles = None
		else:
			self.fam_called = True
			self.fam_alleles = convert_freq_alleles(split_line[7])
	
		if self.rat_cpg == '.':
			self.rat_called = False
			self.rat_alleles = None

		else:
			self.rat_called = True
			self.rat_alleles = convert_freq_alleles(split_line[9])

	def HWE_fail(self,threshold = 0.0002):
	# If a line fails the HWE test, it returns a True
		if self.HWE =='.':
			return False
		elif float(self.HWE) >= threshold:
			return False
		elif float(self.HWE) < threshold:
			return True
			
	def biallelic(self):
		if self.rat_alleles:
			rat = self.rat_alleles
		else:
			rat = self.cast_alleles
		if self.fam_alleles:
			fam = self.fam_alleles
		else:
			fam = self.cast_alleles
		bases = set([i[0] for i in self.cast_alleles] + [i[0] for i in rat] + [i[0] for i in fam] )
		if len(bases) > 2:
			return False
		else:
			return True
	def MAF(self):
		maf = min([i[1] for i in self.cast_alleles])
		if maf == 20:
			return 0
		else:
			return maf
	def diverged(self, outgroup):
		if outgroup not in ['rat','fam']:
			print 'WRONG OUTGROUP GIVEN IN diverged() function'
			return
		if outgroup == 'rat':
			outG = self.rat_alleles
		else:
			outG = self.fam_alleles
		if random_allele(self.cast_alleles)!= random_allele(outG):
			return True
		else: return False


def elementsByGene(fileInput):
	start = True
	for i in gzip.open(fileInput):
		z=i.strip().split()
		
		if start == True:
			currentGeneName = z[12]
			currentGene = [z]
			start = False
			continue
		if z[12] != currentGeneName:
			print currentGeneName
			yield currentGene
			currentGene = [z]
			currentGeneName = z[12]
		elif z[12] == currentGeneName:
			currentGene.append(z)

#def sfsData(intron,freq, depth = 0):
#	chrom = intron[0]
#	start = intron[3]
#	end = intron[4]
#	all_freqs = []
#	nonCpG_freqs = []
#	divergence = Counter()
#
#	for i in freq.fetch(chrom,int(start),int(end)):
#		z = i.strip().split()
#		x = freq_line(z)
#		if x.fam_called and x.cast_called and x.rat_called:pass
#		else:continue
#		if x.HWE_fail():continue
#		if x.depth < depth:continue
#		if not x.biallelic: continue
#		if '1' not in [x.fam_cpg, x.cast_cpg, x.rat_cpg]:
#			nonCpG_freqs.append(x.MAF())
#			if x.diverged('rat'):
#				divergence['ncpgRAT']+=1
#			if x.diverged('fam'):
#				divergence['ncpgFAM']+=1
#		all_freqs.append(x.MAF())
#		if x.diverged('rat'):
#			divergence['allRAT']+=1
#		if x.diverged('fam'):
#			divergence['allFAM']+=1
#
#	all_SFS = SFS_tools.SFS_from_all_frequencies(all_freqs,20)
#	ncpg_SFS = SFS_tools.SFS_from_all_frequencies(nonCpG_freqs,20)
#
#	return {'all_SFS':all_SFS,
#		'ncpg_SFS':ncpg_SFS,
#		'divergence':divergence}	


