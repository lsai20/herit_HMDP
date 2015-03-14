#!/usr/bin/python

# HMDP heritability
# LG
# samples random phenotype from each strain in .pheno file, or uses first phenotype
# assumes all indiv in strain are together in subsequent rows
# handles missing values of -9, any capitalization of nan in phenotype (omitted from average)
# phenotypes where all indiv in phenotype are missing are given value of -9

# input format:
# strain1_0	strain1	pheno0_value	pheno1_value .... phenoN_value

import os
import sys
import random


def usage():
	print("Usage: python strainFirstPheno.py phenotypeFile outputFile numPhenotypes")


def isNonMissingVal(s):
	'''given string s representing a float, return false iff s is a missing value'''
	if s == '-9' or s.lower() == 'nan':
		return False
	return True

if len(sys.argv) == 4:
	fileName = sys.argv[1]
	outputName = sys.argv[2]
	numPheno = int(sys.argv[3])

else: # wrong number of options
    usage()
    sys.exit()


# check if output file already exists to avoid overwriting
if os.path.exists(outputName): 
	print("ERROR: Output file \'%s\' already exists. Please specify a unique name for output." % outputName)
	sys.exit()

isDummyStrain = True
currentStrain = ""
runningTotals = [0.0] * numPheno # for each pheno, sum of phenotype values
numIndivs = [0] * len(runningTotals) # for each pheno, number of indiv without missing value


outf = open(outputName,'w')

with open(fileName,'r') as inf:

	for line in inf:
		lineL = line.strip().split()

		# parse new indiv
		newStrain = lineL[1] # strain is in 2nd col
		newVals = [floatOrMiss(val) for val in lineL[2:]] # 0 if missing, whatever it is otherwise
		newIndivs = [isNonMissingVal(val) for val in lineL[2:]] # True/1 if valid, False/0 if missing

		# if we've reached end of current strain
		if newStrain != currentStrain:
			if isDummyStrain == False:
				# compute avg of old strain and write
				avgPhenos = [str(divideMaybe(x[0],x[1])) for x in zip(runningTotals, numIndivs)]
				avgLine = currentStrain + "_avg\t" + currentStrain + "\t" +  "\t".join(avgPhenos) + "\n"
				outf.write(avgLine)
			
			else: # else just finished dummy strain
				isDummyStrain = False

			# reset vals
			currentStrain = newStrain
			runningTotals = [0.0] * numPheno
			numIndivs = [0] * len(runningTotals)
			

		# update with new individual
		runningTotals = [sum(x) for x in zip(runningTotals, newVals)]
		numIndivs = [sum(x) for x in zip(numIndivs, newIndivs)]


# compute last strain avg after EOF
avgPhenos = [str(divideMaybe(x[0],x[1])) for x in zip(runningTotals, numIndivs)] 
avgLine = currentStrain + "_avg\t" + currentStrain + "\t" +  "\t".join(avgPhenos) + "\n"
outf.write(avgLine)
outf.close()

