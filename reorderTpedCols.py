#!/usr/bin/python

# HMDP heritability
# LG
# reorder indivs in .tped according to order in a phenotype file (where indiv are in rows in correct order)

def usage():
	print("Usage: python reorderTpedCols.py genotypeFileToReorder phenotypeFile outputFile")


if len(sys.argv) != 4:
    usage()
    sys.exit()

genoToReorderFile = sys.argv[1] # .tped with columns to be reordered
phenoRefFile = sys.argv[2] 		# reference for how genotypes will be reordered
outputFile = sys.argv[3]		# output is tped format with order that matches phenoRefFile

orderD = {} 	# orderD[newPos] = oldPos
newPos = -1     # position in ref file will be new pos for corresponding geno
orderL = []		# orderL[newPos] = oldPos

'''
# ref file format
B_0    B    animalID    stuff
B_1    B    animalID    stuff
B_2    B    animalID    stuff


# tped-to-reorder format
# stuff1 stuff2 stuff3 stuff4     geno_old18 geno_old20 
'''

# find reordering
with open(phenoRefFile, 'r') as reff:
	for line in reff:
		lineL = line.split()
		oldPos = int(lineL[2]) # Animal_ID, float cast to int
		newPos += 1
		orderD[newPos] = oldPos

# reorder columns in tped (may be long)
with open(genoToReorderFile) as inf:
	for line in inf:
		oldLineL = line.strip().split()
		newLineL [ for oldPos in orderL] 



