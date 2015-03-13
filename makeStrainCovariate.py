# heritability in HMDP
# LG
# Generate strain matrix from .pheno 
#       (or any file with indiv as row, 2nd col is strain)

import sys


def usage():
    print("Useage: python makeStrainCovariate.py strainPhenoFile outputFileName")

if len(sys.argv) != 3:
    usage()
    sys.exit()


strainPhenoFile = sys.argv[1]   # 2nd col should be strain
outputName = sys.argv[2]

strainStarts = [] # row number of first member in each strain

currentStrain = ''
rowNum = -1 # current row in file

outf = open(outputName, 'w')

with open(strainPhenoFile,'r') as phenoinf:
    for line in phenoinf:
        rowNum += 1
        newStrain = (line.split(None, 3))[1] # second col is strain
        if newStrain != currentStrain:
            strainStarts.append(rowNum)
            currentStrain = newStrain

numIndivs = rowNum + 1
strainStarts.append(numIndivs) # also add index of one past last indiv
print(strainStarts)

with open(outputName, 'r'):

    numStrains = len(strainStarts) - 1
    print(numStrains)
    for i in range(numStrains): # create block of 1's for each strain
        # start positions of current strain and next strain
        start = strainStarts[i]
        startNext = strainStarts[i+1]

        # one row for each member in strain
        for j in range(startNext-start):
            # with one col of 1 for each member
            L = ['0']*(start) + ['1']*(startNext-start) + ['0']*(numIndivs - startNext)
            outf.write(" ".join(L) + "\n")
        


