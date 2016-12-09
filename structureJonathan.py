import numpy as np





# read in the snp data
def read_file():
	individuals = []
	snps = []
	with open("hapmap3.ped") as f:
		for line in f:

			# split on spaces, separate into list of ID's and snp data
			splitLine = line.split()
			individuals.append(splitLine[0:1])
			snps.append(splitLine[6:])

	return individuals, snps


def mcmc_no_admixture(individuals, snps, K):
	number_individuals = len(individuals)

	num_alleles = 2

	# get the number of loci, should be an even number
	num_loci = len(snps[0]) / 2

	# initialize z0 for all individuals from uniform distribution
	z = [np.random.randint(1,K+1) for _ in range(number_individuals)]

	#------ itterate m times  -------- #

	#--------Step 1----------- #

	# set up nklj, a 3D matrix of counts that is indexed by k-> l-> j
	n = np.zeros((K,num_loci,num_alleles))
	# loop through all individuals
	for idx, individual in enumerate(individuals):

		# all the lists are indexed by a common index
		# this particular individual belongs to population k
		k = z[idx]

		# itterate over each locus 
		# loci within snps for each individuals also allign by index

		# just the loci for the current individual
		loci = snps[idx]
		for locus in range(0,num_loci,2):

			# grab allele counts
			allele1 = int(loci[locus])
			allele2 = int(loci[locus+1])

			#print "allele1: {0}".format(allele1)
			#print "allele2: {0}".format(allele2)
			

			# don't count negative counts
			if (allele1 >=0):
				n[k-1][locus][0] = n[k-1][locus][0] + allele1

			# second allele
			if (allele2 >=0):
				n[k-1][locus][1] += allele2

	# fill out matrix p , which stores pkl0, where pkl1 = 1 - pkl0
	p = np.zeros((K,num_loci))
	
	for k in range(0,K):
		for l in range(0,num_loci):

			# sample from dirichlet distribution, return the probablility of allele1
			# store 
			p[k][l] = np.random.dirichlet([1 + n[k][l][0], 1 + n[k][l][1] ])[0]

	print p



def main():
	# read in data
	ids, snps = read_file()
	mcmc_no_admixture(ids,snps,K=3)





if __name__ == "__main__":
	main()