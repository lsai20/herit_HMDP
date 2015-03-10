#!/usr/bin/python

# HMDP heritability
# LG
# finds average phenotype value for each strain from .tfam file
# just a copy of strainAvgPheno.py except the phenotype is in the fifth column instead of starting in the third (PLINK format says it should be in sixth column)
# and assumes one phenotype
# TODO NOT YET WORKING, PROBABLY NOT NEEDED

# assumes all indiv in strain are together in subsequent rows
# handles missing values of -9, any capitalization of nan in phenotype (omitted from average)
# phenotypes where all indiv in phenotype are missing are given value of -9

# input format:
# PLINK tfam: strain1_0	strain1	stuff stuff stuff pheno0_value
# actual input format: strainA_0 strainA 0 2 5.1234 NULL

import sys
import random


def usage():
	print("Usage: python strainAvgTfam.py tfamFile outputFile numPhenotypes")

def isNonMissingVal(s):
	'''given string s representing a float, return false iff s is a missing value'''
	if s == '-9' or s.lower() == 'nan':
		return False
	return True

def floatOrMiss(s):
	'''given string s representing a float, return 0 if s is a missing value. 
	else return value of s'''
	if s == '-9' or s.lower() == 'nan':
		return 0.0
	return float(s)

def divideMaybe(a,b):
	'''divide a/b, return -9 if b is zero'''
	if b == 0:
		return -9
	return a/b


if len(sys.argv) != 4:
    usage()
    sys.exit()

fileName = sys.argv[1]
outputName = sys.argv[2]
numPheno = int(sys.argv[3])
isDummyStrain = True
currentStrain = ""
runningTotals = [0.0] # for each pheno, sum of phenotype values
numIndivs = [0] # for each pheno, number of indiv without missing value

### read first line and count phenotypes
#lineL = inf.readline().strip().split()
#numPheno = len(lineL[2:])
#print("Number of phenotypes found: %d" % numPheno)

outf = open(outputName,'w')

with open(fileName,'r') as inf:

	for line in inf:
		lineL = line.strip().split()

		# parse new indiv
		newStrain = lineL[1] # strain is in 2nd col
		# Change from strainAvgPheno: pheno is in 5th col instead of 3rd, so [4:] instead of [2:]
		newVals = [floatOrMiss(val) for val in lineL[4]] # 0 if missing, whatever it is otherwise
		newIndivs = [isNonMissingVal(val) for val in lineL[4]] # True/1 if valid, False/0 if missing

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

