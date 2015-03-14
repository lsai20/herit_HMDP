#!/usr/bin/python

# HMDP heritability
# LG
# get first indiv in each strain
# assumes all indiv in strain are together in subsequent rows

# input format:
# strain1_0	strain1	pheno0_value	pheno1_value .... phenoN_value

import os
import sys
import random


def usage():
	print("Usage: python strainFirstPheno.py phenotypeFile outputFile")

if len(sys.argv) != 3:
    usage()
    sys.exit()

fileName = sys.argv[1]
outputName = sys.argv[2]

# check if output file already exists to avoid overwriting
if os.path.exists(outputName): 
	print("ERROR: Output file \'%s\' already exists. Please specify a unique name for output." % outputName)
	sys.exit()


def isNonMissingVal(s):
	'''given string s representing a float, return false iff s is a missing value'''
	if s == '-9' or s.lower() == 'nan':
		return False
	return True

def replaceNan(s):
	if s.lower() == 'nan':
		return '-9'
	else:
		return s

outf = open(outputName,'w')
currentStrain = ""

with open(fileName,'r') as inf:

	for line in inf:
		lineL = line.strip().split()
		newStrain = lineL[1] # strain is in 2nd col
		if newStrain != currentStrain:
			# outf.write(line) # just write the unedited line for first indiv

			# replace any nan with -9
			outLineL = [replaceNan(s) for s in lineL]
			outf.write(" ".join(outLineL) + "\n")
			currentStrain = newStrain

outf.close()
