# somaticSmallVariantFilter
Original scripts for Bertha Cancer Workflow 2.0 components "Somatic Variant Filtering" and "Annotate Somatic Small Variant VCF"

## Overview

The first component "Somatic Variant Filtering" calculates, and annotates the INFO column of the somatic small variant VCF with, a Fisher's test phred score for each PASS SNV. This score represents the likelihood the called variant is a false positive due to systematic sequencing/mapping errors. In the cancer workflow the component immediately procedes cellbase annotation.  The second component "Annotate Somatic Small Variant VCF" immediately follows cellbase annotation and adds flags to the vcf FILTER column and annotates the INFO column using the cellbase annotation JSON.  

## Somatic Variant Filtering

Is run on LSF using somaticFisherPhred.sh:
bsub -n 12 -J "noiseFlag" -q bio -P Analysis -o logs/somaticFisherPhred_%J.log -e logs/somaticFisherPhred_%J.err  "bash somaticFisherPhred.sh -b tumourBamFile -v somaticVcf -o outputDirectory/"
NB forward slash required at end of output directory

tumourBamFile: The tumour BAM file generated for running the CANVAS component  
somaticVcf: The normalised small variant vcf
outputDirectory:  location of all output

The final output (annotated vcf) is

outputDirectory/fisherVCF/sampleID.vcf.gz

## Annotate Somatic Small Variant VCF   
 
Is run on LSF using from_annotation_to_vcf.sh:
bsub -n 1 -J "annot" -q bio -P Analysis -o logs/job_%J.log -e logs/job_%J.err "bash from_annotation_to_vcf.sh -v fisherVCF -j annotationJSON -o outputDir/"
NB forward slash required at end of output directory

fisherVCF: vcf generated from component "Somatic Variant Filtering" 
annotationJSON: annotation JSON generated from cellbase 
outputDir: location of all output

The final output (fully annotated vcf) is

outputDirectory/annotatedVCF/sampleID.vcf.gz

## Annotate Somatic Small Variant VCF multiSample

Identical to component AnnotateSomaticSmallVariantVCF but with a small change so that it works with multi sample somatic vcf

NB The componet somaticVariantFiltering works for both single and multi sample somatic vcf without modification.
