import re
import sys

####Add FISHER phred score to vcf####

#Open file containing SNV Fisher Phred score

fileNameFISHER = sys.argv[2]
fisherFile = open(fileNameFISHER, "r")

#Create function to read one line at a time
def read_vars(FISHER_FILE):
    for lineFisher in FISHER_FILE:
        varFisher = lineFisher.rstrip('\n')
        varFisher = re.split(r'\t+',varFisher)
        return varFisher

#Read in first Fisher result
functionRe = read_vars(fisherFile)
if functionRe is None:
    curVar = ["finished", "finished", "finished", "finished", "finished"]
else:
    curVar = functionRe

#Loop through vcf file and add Fisher Phred score to SNV

fileNameVCF = sys.argv[1]

with open(fileNameVCF) as VCF:
    for lineVCF in VCF:
        varVCF = lineVCF.rstrip('\n')

#Print header without modifying
        if varVCF.startswith("#"):
            print(lineVCF.rstrip())

        else:   
            varVCF = re.split(r'\t+',varVCF)
#If not SNV then print
            if varVCF[3] not in ["A", "C", "G", "T"] or varVCF[4] not in ["A", "C", "G", "T"]:
                print("\t".join(varVCF))
#Check if chromosome, position, ref and alt match, and if they do add the Fisher's test score
            elif (curVar[0] == varVCF[0] and curVar[1] == varVCF[1] and curVar[3] == varVCF[3] and curVar[4] == varVCF[4]):
                varVCF[7] = varVCF[7] + ";SomaticFisherPhred=" + str(curVar[7]) 
                print("\t".join(varVCF))
                
                functionRe = read_vars(fisherFile)
                if functionRe is None:
                    curVar = ["finished", "finished", "finished", "finished", "finished"]
                else:
                    curVar = functionRe                

            else:
                print("\t".join(varVCF))        
