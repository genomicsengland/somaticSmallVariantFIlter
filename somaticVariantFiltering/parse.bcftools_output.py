import re
import sys

####Add allele depth pileup to vcf####

#Open bcftools pileup (BTP) file and create function to read by line 

fileNameBTP = sys.argv[2]
BTPfile = open(fileNameBTP, "r")

def read_vars(BTP_FILE):
    for lineBTP in BTP_FILE:
        varBTP = lineBTP.rstrip('\n')
        varBTP = re.split(r'\t+',varBTP)
        return varBTP

curVar = read_vars(BTPfile)

#Loop through vcf file and add allele depth

fileNameVCF = sys.argv[1]

with open(fileNameVCF) as VCF:
    for lineVCF in VCF:
        varVCF = lineVCF.rstrip('\n')
        varVCF = re.split(r'\t+',varVCF)


#Check if chromosome and position match
#If yes add alelle depths  

        if (curVar[0] == varVCF[0] and curVar[1] == varVCF[1]):
            alleles = curVar[2].split(',')
            tumourAD = map(int, curVar[3].split(','))            
            try:
                ref = alleles.index(varVCF[3])
                tumourRefD = tumourAD[ref+1]
            except ValueError:
                tumourRefD = 0
            try:
                alt = alleles.index(varVCF[4])
                tumourAltD = tumourAD[alt+1]
            except ValueError:
                tumourAltD = 0
            tumourTotAD = sum(tumourAD)
            tumourADout = "%s,%s:%s" % (tumourRefD,tumourAltD,tumourTotAD)
            varOut=varVCF[:8]
            varOut.extend(("AD:DP",tumourADout))
            print("\t".join(varOut))


#Read in next pileup position and check if end of file

            functionRe = read_vars(BTPfile)
            if functionRe is None:
                curVar = ["finished", "finished"]
            else:
                curVar = functionRe

#If not then allele depths are zero

        else:
            varOut=varVCF[:8]
            varOut.extend(("AD:DP","0,0:0"))   
            print("\t".join(varOut))
