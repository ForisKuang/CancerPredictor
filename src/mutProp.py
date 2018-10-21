import sys, csv

def main():
	args = sys.argv
	if (len(args) < 3):
		print("Usage: ./mutProp.py <total gene lengths file> <indices file> <OPTIONAL output filename>")
		exit(-1)

	if (len(args) < 4):
		outFileName = "geneProps.csv"
	else:
		outFileName = args[3]


	geneFile = open(args[1], 'r')
	genesLen = list(csv.reader(geneFile))
	# Figure out relavent genes and total lengths
	geneLens = {}
	for i in range(1, len(genesLen)):
		line = genesLen[i]
		geneLens[line[0]] = int(line[2])

	genes = geneLens.keys()

	indsFile = open(args[2], 'r')
	indsLen = list(csv.reader(indsFile))
	# Count number of relavent rows with genes
	indLens = {}
	for i in range(1, len(indsLen)):
		line = indsLen[i]
		gene = line[0]
		if (gene in geneLens):
			if (not gene in indLens):
				indLens[gene] = 0

			if (int(line[2]) > 0):
				indLens[gene] += 1

	outFile = open(outFileName, 'w+')
	for gene in genes:
		output = gene + ',' + str(100*indLens[gene]/geneLens[gene]) + "%\n"
		outFile.write(output)


if __name__ == "__main__":
	main()