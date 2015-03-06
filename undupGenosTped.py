#!/usr/bin/python


# heritability in HMDP
# LG
# unduplicate genotypes while preserving order of strains in phenotype file


import sys

sys.argv


def usage():
    print("Useage: python undupGenosTped.py strainPhenoFile dupGenoFile outputFileName")

if len(sys.argv) != 4:
    usage()
    sys.exit()


strainPhenoFile = sys.argv[1]   #
dupGenoFile = sys.argv[2]       # duplicated tped
outputName = sys.argv[3]

keepIndivs = [] # row number of first member in each strain

currentStrain = ''
rowNum = -1 # current row in file

# find row id of first indiv in each strain
with open(strainPhenoFile,'r') as phenoinf:
    for line in phenoinf:
        rowCount += 1
        newStrain = (line.split('\s', 3))[2]
        if newStrain != currentStrain:
            keepIndivs.append(rowCount)

# write tped with only one geno column per indiv
with open(dupGenoFile,'r') as genoinf, open(outputName,'w') as outf:
    for line in genoinf:
        snpL = line.strip().split()
        # first four cols are chr/snp id stuff, rest are genotype for each indiv
        undupL = snpL[:4] + [snpL[4:][i] for i in keepIndivs]
        outf.write("\t".join(undupL) + "\n")








