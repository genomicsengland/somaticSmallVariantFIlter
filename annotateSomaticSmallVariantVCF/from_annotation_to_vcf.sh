#!/bin/bash

# Submit with 
# bsub -n 1 -J "annot" -q bio -P Analysis -o logs/job_%J.log -e logs/job_%J.err "bash from_annotation_to_vcf.sh -v fisherVCF -j annotationJSON -o outputDir"
# e.g.
# bsub -n 1 -J "annot" -q bio -P Analysis -o logs/job_%J.log -e logs/job_%J.err "bash from_annotation_to_vcf.sh -v /home/jmitchell1/noiseModelBertha_passArg/annotateSomaticSmallVariantVCF/testOutput/fisherVCF/LP3000396-DNA_G02.vcf.gz -j /home/jmitchell1/noiseModelBertha_passArg/annotateSomaticSmallVariantVCF/testInput/LP3000396-DNA_G02.json.gz -o /home/jmitchell1/noiseModelBertha_passArg/annotateSomaticSmallVariantVCF/testOutput/" 


#Set the cancer test environment to replicate Bertha
source /genomes/software/src/test-venvs/cancer-test/bin/activate

#Load modules
module load bcftools/1.9

#If logs directory doesn't exist then make it
mkdir -p logs

#Arguments are file paths to fisher somatic small variant VCF
#json annotation file
#and output directory
while getopts v:j:o: option
do
case "${option}"
in
v) vcf=${OPTARG};;
j) json=${OPTARG};;
o) outDir=${OPTARG};;
esac
done



mkdir -p ${outDir}
samplename=${vcf##*/}
samplename=${samplename%%.*}


#Create directories and file names
annotated_vcf_dir="${outDir}annotatedVCF/"
mkdir -p ${annotated_vcf_dir}
annotated_vcf="${annotated_vcf_dir}${samplename}.vcf.gz"

#Set running directory to directory of this script
DIR=`dirname $0`

#Create VCF with additional filters and annotation, and additional annotations (python script takes sample name from vcf, everything before first .)
python ${DIR}/from_annotation_to_vcf.py -v ${vcf} -j ${json} -o ${annotated_vcf_dir}
tabix -p vcf ${annotated_vcf}
