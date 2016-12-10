import numpy as np
import functools, operator
from math import log, exp
from numpy import log1p
import sys

# ------- CONSTANTS ----------#
ALLELE_1 = 1
ALLELE_2 = 2




#----------- Helper Functions --------------#

# Return the sum of two log probabilities
# requires two log probabilities a and b
def sumLogProb(a, b):
	if a > b:
		return a + log1p(exp(b - a))
	else:
		return b + log1p(exp(a - b))





# helper function for computing the probability 
def pHelper(p, k):
	return reduce(operator.add, p[k])



# read in the snp data
def read_file(fileName):
	individuals = []
	snps = []
	with open(fileName) as f:
		for line in f:

			# split on spaces, separate into list of ID's and snp data
			splitLine = line.split()
			individuals.append(splitLine[0:1])
			snps.append(map(int,splitLine[6:]))

		# convert all snps to integers
		snps = np.array(snps)
	return individuals, snps


def mcmc_no_admixture(individuals, snps, K):
	number_individuals = len(individuals)
	print "number of individuals: {0}".format(number_individuals)


	num_alleles = 2

	# get the number of loci, should be an even number
	num_loci = len(snps[0]) / 2
	print "number of loci: {0}".format(num_loci)

	# initialize z0 for all individuals from uniform distribution
	z = np.array([np.random.randint(0,K) for _ in range(number_individuals)])
	#z = [2, 1, 2, 1, 2, 2]
	
	'''
	Algorithm overview


	Terminate when likelyhood of data is constant

	'''
	#------ itterate m times  -------- #
	for m in range(10):
		print "itteration {0}".format(m)

		#--------Step 1----------- #
		print "Run Step 1"
		# set up nklj, a 3D matrix of counts that is indexed by k-> l-> j
		n = np.ones((K,num_loci,num_alleles))
		# loop through all individuals
		for idx, individual in enumerate(individuals):

			# all the lists are indexed by a common index
			# this particular individual belongs to population k
			k = z[idx]

			# itterate over each locus 
			# loci within snps for each individuals also allign by index

			# just the loci for the current individual
			loci = snps[idx]

			# itterate over each locus
			for locus in range(0,num_loci):

				# grab chromosomes at this locus
				chromosome1, chromosome2 = int(loci[2*locus]), int(loci[2*locus+1])

				#print "allele1: {0}".format(allele1)
				#print "allele2: {0}".format(allele2)
				
				# count alleles 
				for chromosome in [chromosome1, chromosome2]:
					if (chromosome == ALLELE_1):
						n[k-1][locus][0] += 1
					elif (chromosome == ALLELE_2):
						n[k-1][locus][1] += 1


		# fill out matrix p , which stores pkl0, where pkl1 = 1 - pkl0
		# use pseudocounts, minimum count is 1
		p = np.ones((K,num_loci))
		
		for k in range(0,K):
			for l in range(0,num_loci):

				# sample from dirichlet distribution, return the probablility of allele1
				# store as logs
				p[k][l] = log(np.random.dirichlet([1 + n[k][l][0], 1 + n[k][l][1] ])[0])

		print "population distribution: {0}".format(z)
		print "matrix: n {0}".format(n)
		print "matrix: p  {0}".format(p)
		#--------Step 2----------- #
		'''
		Description of Step 2
		We simulate z(i) from Equation A8. We need to calculate the log likelyhood
		 of all the loci for each individual. This will serve as the normalization
		 factor. an individual will be sampled from a population k by a weighted
		 vector [a(k1}/norm, a(k2)/norm, a(k3)/norm] where a is Pr(x(i)|P, z(i)=k) 
		'''
		print "Run Step 2"

		# itterate over each individual
		for idx, individual in enumerate(individuals):

			# just the loci for the current individual
			loci = snps[idx]
			kDistribution = np.array( probHelper(loci, num_loci, p))
			sum_k_dist = reduce(sumLogProb, kDistribution)
			# print "before normalizing: {0}".format(kDistribution)


			# normalize
			kDistribution = map(exp, kDistribution - sum_k_dist)
			print "after normalizing: {0}".format(kDistribution)



			# sample from this weighted distribution, assign new k to this individual:
			z_new = np.random.choice(K, 1,p=kDistribution)[0]
			if z[idx] != z_new:
				print "individual {0} changed from {1} to {2}".format(idx, z[idx], z_new)
				z[idx] = z_new






# calculate the product of p's based on equation A8, according to the alleles
# at the particular locus
# using log probabilities
def alleleHelper(pLocus, allele1, allele2):

	# terms = [0,0]
	# for idx, allele in enumerate([allelle1, allele2]):
	# 	if allele == 1:
	# 		terms[idx] = pLocus
	# 	elif allele == 2:
	# 		terms[idx] = log(1 - exp(pLocus))
	# 	else:
	# 		terms[idx] = 0
	# return sum(terms)
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



# compute Pr(x(i)|P, z(i)=k) for all k
def probHelper(loci, num_loci, p):
	# return [ sum([ alleleHelper(pLocus,loci[2*idxLocus],loci[2*idxLocus+1])  for idxLocus, pLocus in enumerate(pk)]) for pk in p]
	kDist = []
	for pk in p:
		
		accum = 0
		for idxLocus, pLocus in enumerate(pk):
			accum += alleleHelper(pLocus, loci[2*idxLocus], loci[2*idxLocus+1])

		kDist.append(accum)

	return kDist


# def logLikelyHood(p):
# 	accum = 0 
# 	for k in 



def main():
	# read in data
	ids, snps = read_file("hapmap3.ped")
	#ids, snps = read_file("smallTest.ped")


	mcmc_no_admixture(ids,snps,K=3)





if __name__ == "__main__":
	main()