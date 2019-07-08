#!/usr/bin/env Rscript


# Arguments: tumour pileup vcf; output vcf
args = commandArgs(trailingOnly=TRUE)

# PON files path
# 7K Indv PoN
PON_dir="/home/jmitchell1/somaticMutations_falsePositives/PoN_7k/chromosomes/"

nreadlines=1000

tumour_germline_file=args[[1]]
outputfile=args[[2]]

chrlist=paste0("chr",c(1:22,"X","Y"))

library("stringr")

# Function to caluclate phred score from pValue

phred=function(x){if(x==0){10000} else {-10*log10(x)}}

# Function to retrieve allele depths from PoN

ExtractPON=function(chr,pos,altbase){
  cmd=paste0("tabix ",PON_dir,chr,"_pos_ref_alt_AD.txt.gz ",chr,":",pos,"-",pos," | cut -f 5")
  stm=system(cmd,intern = T)
  if(length(stm)==0){return("VariantNotFound")}
names(stm)=c("A","C","G","T")
if(stm[[altbase]]=="0,0:0,0:0,0"){
  mx=max(as.integer(str_split_fixed(stm,",|:",3)[,2]))
  return(paste0("0,",mx,":0,",mx,":0,0"))
}
return(stm[[altbase]])
}

# Start of ProcessVCFlines function

ProcessVCFlines=function(lines){
vcf=str_split_fixed(lines,"\t",10)
print(unique(vcf[,1]))

vcf=vcf[vcf[,1] %in% chrlist,,drop=F]

nlines=nrow(vcf)

if(nlines==0){return(vcf)}

rc=array(as.integer(str_split_fixed(vcf[,10],",|:",3)),dim=c(nrow(vcf),3))

for (j in 1:nlines){
  
  pon=ExtractPON(vcf[j,1],vcf[j,2],vcf[j,5])
  if(pon=="VariantNotFound"){
  cnt=c(0,0,0,0,0,0)
  } else {
  cnt=as.integer(str_split_fixed(pon,",|:",6))
  }
 
# Perform Fisher test
  x=matrix(c(rc[j,2],rc[j,3]-rc[j,2],cnt[3],cnt[4]-cnt[3]),nrow=2)
  f=fisher.test(x,alternative="greater")
  fisher_00=phred(f$p.value)  
  vcf[j,8]=round(fisher_00)
}
  return(vcf)

}

# End of ProcessVCFlines function

# Processing

inputpipe=file(tumour_germline_file,"rt",blocking=T)
if(file.exists(outputfile)){file.remove(outputfile)}

print("STARTED!")
print(Sys.time())
j=0
while (TRUE){
  j=j+1
  lines=readLines(inputpipe,nreadlines)
  if(length(lines)==0){break}
  print(paste("Processing the ",j," set of ",nreadlines," lines"))
  print(Sys.time())
  res=ProcessVCFlines(lines)
  if(nrow(res)>0){
  write.table(res,
              file = outputfile,
              append = T,
              col.names = F,row.names = F,quote = F,sep="\t")
  }
  else {print("No vars collected")}
  
}

close(inputpipe)
print("FINISHED!")
print(Sys.time())
