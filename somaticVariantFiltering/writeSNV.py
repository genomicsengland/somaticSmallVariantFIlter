import re
import gzip
import sys

####Pull out auto/sex SNVs which are PASS####

#Input VCF 
fileNameVcf = sys.argv[1]
#Output file
outFileName = sys.argv[2]
snvVars = open(outFileName,'w')

autosomeSex = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
bases = ["A", "C", "G", "T"]

with gzip.open(fileNameVcf,'r') as vcfFile:
    for linevcf in vcfFile:
        varVCF = linevcf.rstrip('\n')

#Skip header
        if varVCF.startswith("#"):
            continue

        else:   
            varVCF = re.split(r'\t+',varVCF)

#If variant is PASS, an SNV and in autosome/sex chromosome then print
            if varVCF[6] == "PASS" and varVCF[0] in autosomeSex and varVCF[3] in bases and varVCF[4] in bases:
                snvVars.write("\t".join(varVCF))
                snvVars.write("\n")
