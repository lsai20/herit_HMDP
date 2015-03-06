#!/usr/bin/python


# heritability in HMDP
# LG
# unduplicate genotypes while preserving order of strains in phenotype file
# really slow, don't use
 
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

print("Now finding individuals to keep in genotype file...")
# find row id of first indiv in each strain
with open(strainPhenoFile,'r') as phenoinf:
    for line in phenoinf:
        rowNum += 1
        newStrain = (line.split(None, 3))[1] # second col is strain
        if newStrain != currentStrain:
            keepIndivs.append(rowNum)

print("Now trimming genotype file...")
snpcount = 0 # print progress
# write tped with only one geno column per indiv
outf = open(outputName,'w') # py2.6 can only do one file in with i believe
with open(dupGenoFile,'r') as genoinf:
    for line in genoinf:
        snpcount += 1
        if snpcount % 2000 == 0:
            print("writing snp %d..." % snpcount)
        snpL = line.strip().split()
        # first four cols are chr/snp id stuff, rest are genotype for indiv
        undupL = snpL[:4] + [snpL[4:][i] for i in keepIndivs]
        outf.write("\t".join(undupL) + "\n")

outf.close()
print("Done unduplicating")






