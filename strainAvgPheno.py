# HMDP heritability
# LG
# finds average phenotype value for each strain from .pheno file
# assumes all indiv in strain are together in subsequent rows
# input format:
# strain1_0	strain1	pheno0_value	pheno1_value .... phenoN_value

# TODO remove nan and -9. sum of float('nan') and any float is nan, and -9 also need to be treated as missing values.
# that individual also needs to be treated as not incrementing the count.

import sys
import random

def usage():
	print("Usage: python strainAvgPheno.py phenotypeFile outputFile")

def sumLines(floatL, line):
	'''adds values in line to floatL, ignoring nan and -9'''


if len(sys.argv) != 3:
    usage()
    sys.exit()

fileName = sys.argv[1]
outputName = sys.argv[2]


with open(fileName,'r') as inf, open(outputName,'w') as outf:
	# read first line
	lineL = inf.readline().strip().split()
	print("first line")
	print(lineL)
	numIndiv = 1 # number of indiv in strain so far
	currentStrain = lineL[1]
	runningTotals = [float(val) for val in lineL[2:]] # totals for each phenotype
	print("Number of phenotypes found: %d" % len(runningTotals))
	
	# rest of indivs
	for line in inf:
		lineL = line.strip().split()
		newStrain = lineL[1] # second column is strain name
		newPhenoVals = [float(val) for val in lineL[2:]]

		if newStrain != currentStrain:
			# find avg of old strain
			print(currentStrain)
			avgPhenos = [str(1.0*val/numIndiv) for val in runningTotals]
			indivID = currentStrain + "_avg"
			avgLine = "\t".join([indivID, currentStrain] + avgPhenos)
			outf.write(avgLine + "\n")

			# reset values before continuing
			runningTotals = [0 for val in runningTotals]
			numIndiv = 0
			currentStrain = newStrain

		numIndiv += 1
		runningTotals = [sum(x) for x in zip(runningTotals, newPhenoVals)]


# TODO compute last strain avg after EOF
avgPhenos = [1.0*val/numIndiv for val in runningTotals]
indivID = currentStrain + "_avg"
avgLine = "\t".join([indivID, currentStrain] + avgPhenos)
open(outfName, 'a').write(avgLine + "\n")

