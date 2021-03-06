import numpy as np
import functools, operator
from math import log, exp
from numpy import log1p
import sys

# ------- CONSTANTS ----------#
ALLELE_1 = 1
ALLELE_2 = 2
NUM_ALLELES = 2




#----------- Helper Functions --------------#

# Return the sum of two log probabilities
# requires two log probabilities a and b
def sumLogProb(a, b):
	if a > b:
		return a + log1p(exp(b - a))
	else:
		return b + log1p(exp(a - b))

'''
read in the snp data
singleLine is a parameter to read in SNP data that's formated by 1 individual / line
Hapmap3 data is in the following format:
1 individual per line
Locusi_1 Locusi_2 ....  

'''
def read_file(fileName, singleLine=True):
	individuals = []
	snps = []
	with open(fileName) as f:
		for line in f:

			# split on spaces, separate into list of ID's and snp data
			splitLine = line.split()
			individuals.append(splitLine[0:1])
			snps.append(map(int,splitLine[6:]))

		# convert all snps to integers
		#snps = np.array(snps)
		number_individuals = len(individuals)
		print "number of individuals read: {0}".format(number_individuals)
		print "number of loci read: {0}".format(len(snps[0])/2)

	return individuals, snps


# --BEGIN no admixture helpers #

''' 
alculate the product of p's based on equation A8, according to the alleles
at the particular locus
requires a log probability pLocus and two allele copies allele1 and allele2
'''
def alleleHelper(pLocus, allele1, allele2):

	if allele1 == 1:
		term1 = pLocus
	elif allele1 == 2:
		term1 = log(1 - exp(pLocus))
	else:
		term1 = 0

	if allele2 == 1:
		term2 = pLocus
	elif allele2 == 2:
		term2 = log(1 - exp(pLocus))
	else:
		term2 = 0

	return term1 + term2
'''
compute Pr(x(i)|P, z(i)=k) for all k
requires a list loci of all of the loci for an individual i
'''
def probHelper(loci, num_loci, p):

	kDist = []
	for pk in p:
		
		accum = 0
		for idxLocus, pLocus in enumerate(pk):
			accum += alleleHelper(pLocus, loci[2*idxLocus], loci[2*idxLocus+1])

		kDist.append(accum)

	return kDist



# --END no admixture helpers #

# --- Admixture Helpers ---#
# map genotype observed at allele to correct probability
# p only stores probability of ALLELE_1
# requires a log probability and an allele of value either ALLELE_1, ALLELE_2, or 0, indicating missing
def single_allele_helper(logProb, allele):
	if allele == ALLELE_1:
		return logProb
	elif allele == ALLELE_2:
		return log(1 - exp(logProb))

	# missing data case
	# not 100% sure what to do here
	else:
		return 0
#-----------END Helper Functions --------------#

'''
STRUCTURE no Admixture
requires a list of individuals, where each element of individuals is a 2 element list 
requires a list of snps, where each element is a list of allele copies of length number_loci*2
(2 allele copies per locus)
'''
def mcmc_no_admixture(individuals, snps, K):
	
	number_individuals = len(individuals)


	# get the number of loci, should be an even number
	num_loci = len(snps[0]) / 2

	print "number of individuals for alg: {0}".format(number_individuals)
	print "number of loci for alg: {0}".format(num_loci)

	# initialize z0 for all individuals from uniform distribution
	z = np.array([np.random.randint(0,K) for _ in xrange(number_individuals)])
	#z = [2, 1, 2, 1, 2, 2]
	
	'''
	Algorithm overview
	In no admixture case, assume all individuals come from 1 population 
	Initialize z with uniform distirbution(K)
	Step 1: Generate matrix nklj, sample p[k][l] from Dirichlet
	Step 2: sample z[i] from A8

	'''

	'''
	Data structure overview (no admixture)
	List indeces map to individuals
	
	z -> list of individuals' k assignment; 1 k per individual. length(z) = # ind.
	snps[i][a] -> a 2D list indexed by individual and allele copy
		length(snps) = number of individuals
		length(snps[i]) for all = # loci*2, since 2 allele copies per locus

	n[k][l][j] -> 3D list of counts of allele j at locus l observed in pop k

	p[k][l] -> 2D list. stores the log probability of observing ALLELE_1. Get probability
	of ALLELE_2 by taking log(1 - exp(p[k][l])) 

	'''

	#------ itterate m times  -------- #
	for m in xrange(50):

		#--------BEGIN Step 1----------- #
		# print "Run Step 1"
		# set up nklj, a 3D matrix of counts that is indexed by k-> l-> j
		n = np.ones((K,num_loci,NUM_ALLELES))

		# loop through all individuals
		for idx, individual in enumerate(individuals):

			# all the lists are indexed by a common index
			# this particular individual belongs to population k
			k = z[idx]


			# loci within snps for each individuals also allign by index

			# just the loci for the current individual
			loci = snps[idx]

			# itterate over each locus
			for locus in xrange(0,num_loci):

				# grab allele copies at this locus
				a_copy_1, a_copy_2 = int(loci[2*locus]), int(loci[2*locus+1])

				# count allele copies 
				for a_copy in [a_copy_1, a_copy_2]:
					if (a_copy == ALLELE_1):
						n[k][locus][0] += 1
					elif (a_copy == ALLELE_2):
						n[k][locus][1] += 1


		# fill out matrix p , which stores pkl0, where pkl1 = 1 - pkl0
		p = np.zeros((K,num_loci))
		
		for k in xrange(0,K):
			for l in xrange(0,num_loci):

				# sample from dirichlet distribution, return the probablility of allele1
				# store as logs
				p[k][l] = log(np.random.dirichlet([1 + n[k][l][0], 1 + n[k][l][1] ])[0])

		print "itter: ({0}) p log likelihood (simpl)  {1}".format(m, p.sum())
		#--------END Step 1----------- #

		#--------BEGIN Step 2----------- #
		'''
		Description of Step 2
		We simulate z(i) from Equation A8. We need to calculate the log likelyhood
		 of all the loci for each individual. This will serve as the normalization
		 factor. an individual will be sampled from a population k by a weighted
		 vector [a(k1}/norm, a(k2)/norm, a(k3)/norm] where a is Pr(x(i)|P, z(i)=k) 
		'''
		individuals_changed = 0

		# itterate over each individual
		for idx, individual in enumerate(individuals):

			# just the loci for the current individual
			loci = snps[idx]
			kDistribution = np.array( probHelper(loci, num_loci, p))
			sum_k_dist = reduce(sumLogProb, kDistribution)
			# print "before normalizing: {0}".format(kDistribution)

			# normalize
			kDistribution = map(exp, kDistribution - sum_k_dist)
			#print "after normalizing: {0}".format(kDistribution)

			# sample from this weighted distribution, assign new k to this individual:
			z_new = np.random.choice(K, 1,p=kDistribution)[0]
			if z[idx] != z_new:
				#print "individual {0} changed from {1} to {2}".format(idx, z[idx], z_new)
				individuals_changed += 1
				z[idx] = z_new

		# print "number of individuals reassigned in itteration {0}: {1}".format(m, individuals_changed)
	print "k=0 ind: {0}".format(np.count_nonzero(z == 0))
	print "k=1 ind: {0}".format(np.count_nonzero(z == 1))
	print "k=2 ind: {0}".format(np.count_nonzero(z == 2))
	#print "final ancestry assignment: {0}".format(z)

# STRUCTURE with admixture
# details described in report
# same input as mcmc_noadmixture
def mcmc_admixture(individuals, snps, K=3):
	
	# get the number of loci, should be an even number
	num_loci = len(snps[0]) / 2

	# get number of individuals
	number_individuals = len(individuals)

	print "number of individuals for alg: {0}".format(number_individuals)
	print "number of loci for alg: {0}".format(num_loci)

	# each observed allele copy xl(i,a) originated in some unkown population zl(i,a)
	# initialize z0 for all alleles for all individuals from uniform distribution
	z = np.array([ [np.random.randint(0,K) for _ in xrange(num_loci * 2)] for _ in xrange(number_individuals)])
	

	# itterate itter times
	for itter in xrange(5000):
		print "running itteration {0}".format(itter)

		#--------Step 1----------- #
		# print "Run Step 1"
		# set up nklj, a 3D matrix of counts that is indexed by k-> l-> j
		# pseudocounts for n
		n = np.ones((K,num_loci,NUM_ALLELES))

		# m -> number of allele copies in indvidual i that orginated in population k
		m = np.zeros((number_individuals, K))

		# q -> proportion of admixture
		q = np.zeros((number_individuals, K))

		# itterate over each individual
		for idx, individual in enumerate(individuals):

			# itterate over each locus 
			# loci within snps for each individuals also allign by index

			# just the loci for the current individual
			# loci stores the diploid SNPS for each locus, so len(loci) should be 2*num_loci
			loci = snps[idx]

			# itterate over each locus
			for locus in xrange(0,num_loci):

				# grab allele copies at this locus
				a_copy_1, a_copy2 = int(loci[2*locus]), int(loci[2*locus+1])

				# grab the populations of origin of the allele copies (of a_copy_1 and a_copy2)
				k1, k2 = z[idx][2*locus], z[idx][2*locus+1]

				#print "allele1: {0}".format(allele1)
				#print "allele2: {0}".format(allele2)
				
				# count allele copy 1 
				if (a_copy_1 == ALLELE_1):
					n[k1][locus][0] += 1
				elif (a_copy_1 == ALLELE_2):
					n[k1][locus][1] += 1

				# count allele copy 2
				if (a_copy2 == ALLELE_1):
					n[k2][locus][0] += 1
				elif (a_copy2 == ALLELE_2):
					n[k2][locus][1] += 1

				# keep track of allele copies seen
				m[idx][k1] += 1
				m[idx][k2] += 1


		#------ sample p-----#
		# fill out matrix p , which stores pkl0, where pkl1 = 1 - pkl0
		# use pseudocounts, minimum count is 1 (should I do this?)
		p = np.zeros((K,num_loci))
		
		for k in xrange(0,K):
			for l in xrange(0,num_loci):

				# sample from dirichlet distribution, return the probablility of allele1
				# store as logs
				p[k][l] = log(np.random.dirichlet([1 + n[k][l][0], 1 + n[k][l][1] ])[0])

		#-----sample q------#
		# for each individual, update q
		for i in xrange(number_individuals):
			# sample from Dirichlet distribution using m
			# alpha = 1
			q[i] = np.random.dirichlet([1 + m[i][0], 1 + m[i][1], 1 + m[i][2] ])

		#--------END Step 1----------- #

		#--------BEGIN Step 2----------- #
		# need to choose a k for each allele copy in each locus, for every individual 

		# simulate zl(i,a)

		# itterate over all z entries
		for i in xrange(number_individuals):
			for locus in xrange(num_loci):
				for a in xrange(locus*2,locus*2+1):

					# generate distribution
					dist = np.array([ log(q[i][k]) + single_allele_helper(p[k][locus], snps[i][a]) for k in xrange(0,K)])
					sum_dist = reduce(sumLogProb, dist)
			
					# normalize
					dist = map(exp, dist - sum_dist)

					# sample zl(i,a) from this distribution
					# assign new k to this allele copy
					z[i][a] = np.random.choice(K, 1,p=dist)[0]

		# --- Stats for each run --- #
		print "p log likelihood (simpl)  {0}".format(p.sum())
		if (itter%20 == 0):
			print "admix fractions after {0} runs : {1}".format(itter,q)
	print "admix fractions after {0} runs : {1}".format(itter,q)
	#--------END Step 2----------- #


def main():
	# read in data
	ids, snps = read_file("hapmap3.ped")
	ids = np.array(ids)
	snps = np.array(snps)
	#ids, snps = read_file("thrush", False)
	#ids, snps = read_file("smallTest.ped")

	# run no_admixture many times
	# for run in xrange(0,30):
	# 	print "run number {0}".format(run)
	#mcmc_no_admixture(ids,snps,K=3)


	mcmc_admixture(ids,snps)

if __name__ == "__main__":
	main()