



# compare outputs of STRUCTURE and jonathanStructure
def main():
	jonathan = []
	with open('s_j_3') as f:
		for line in f:
			jonathan.append(map(int,line.split()))

	structure = []
	with open('f4_Q') as f:
		for line in f:
			structure.append(map(float,line.split()))


	jonathanOut = []
	for line in jonathan:
		if line[0] == 1:
			jonathanOut.append("Red")
		elif line[1] == 1:
			jonathanOut.append("Blue")
		elif line[2] == 1:
			jonathanOut.append("Green")
	structureOut = []
	for line in structure:
		if line[0] == 1:
			structureOut.append("Red")
		elif int(line[1]) == 1:
			structureOut.append("Green")
		elif int(line[2]) == 1:
			structureOut.append("Blue")
		else:
			print "somethign went wrong"

	error = 0
	for i in range(len(jonathanOut)):
		if jonathanOut[i] != structureOut[i]:
			error +=1
	print error


if __name__ == "__main__":
	main()