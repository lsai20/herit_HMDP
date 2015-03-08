#!/usr/bin/python

# HMDP heritability
# LG
# NOTE: Not finished. Use undupGenosTped on duplicated instead.
# reorder indivs in .tped according to order in two reference phenotype/tfam file (where indiv are in rows in original order and in correct order)

import sys

def usage():
	print("Usage: python reorderTpedCols.py genotypeFileToReorder oldOrderFile newOrderFile outputFile")


if len(sys.argv) != 5:
    usage()
    sys.exit()

genoToReorderFile = sys.argv[1] # .tped with columns to be reordered
oldOrderFile = sys.argv[2]      # ref for how genotype file is currently ordered
newOrderFile = sys.argv[3] 		# ref for how genotypes will be reordered
outputFile = sys.argv[4]		# output is tped format with order that matches phenoRefFile

orderD = {} 	# orderD[newPos] = oldPos
newPos = -1     # position in ref file will be new pos for corresponding geno
orderL = []		# orderL[newPos] = oldPos

'''
# pheno or tfam file format
stuff1    strainA    stuff stuff stuff ...
stuff1    strainB    stuff ...
stuff1    strainC    stuff ...


# tped-to-reorder format
# stuff1 stuff2 stuff3 stuff4     geno_old18 geno_old20 
'''

# find old strain order
strainD = {}

with open(oldOrderFile, 'r') as origf:
    strainsL_oldOrder = [line.split()[1] for line in origf.readlines()]

with open(oldOrderFile, 'r') as origf:
    pos = -1 # current position
    oldStrainL = []
    for line in origf:
        strain = line.split()[1]
        oldStrainL.append(strain)
        pos += 1
        strain = line.split()[1]
        strainD[strain] = (-1, pos)

# find new strain order
with open(newOrderFile, 'r') as newOrderf:
    strainsL_newOrder = [line.split()[1] for line in newOrderf.readlines()]

with open(newOrderFile, 'r') as newOrderf:
    pos = -1
    newStrainL = []
    for line in newOrderf:
        strain = line.split()[1]
        newStrainL.append(strain)
        pos += 1
        (dummy,oldPos) = strainD[strain]
        strainD[strain] = (pos, oldPos)


# map from to new order to old order
# oldPosOf[newPos] = oldPos
for key, val in strainD.items():
    print(key,val)

print(sorted(oldStrainL))
print(sorted(newStrainL))

oldPosOf = sorted(strainD.values())
#print(oldPosOf)

posMapD = {}
for newPos, oldPos in strainD.values():
    posMapD[newPos] = oldPos

# reorder columns in tped (may be long)
with open(genoToReorderFile) as inf:
    for line in inf:
        lineL = line.strip().split()
        oldGenoL = lineL[4:]
        numGeno = len(oldGenoL)
        newGenoL = [oldGenoL[posMapD[newPos]] for newPos in range(numGeno)]
        sys.stdout.write("\t".join(lineL[:4] + newGenoL) + "\n")


