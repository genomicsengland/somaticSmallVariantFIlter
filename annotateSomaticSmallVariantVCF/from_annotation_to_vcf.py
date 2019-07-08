# Requirements: https://pypi.gel.zone/genomics/dev/pythoncommonlibs/0.14.2 Documentation for Annotation
# Models: http://gelreportmodels.genomicsengland.co.uk/html_schemas/org.opencb.biodata.models.variant.avro/1.3.0-SNAPSHOT/variantAnnotation.html#/schema/org.opencb.biodata.models.variant.avro.VariantAnnotation

import cStringIO
import argparse
import sys
import subprocess
import os
from re import sub
from DataProducers.VariantReaders import JsonGelVariantReader
from DataProducers.VariantProducers import GVCF
import inspect

def add_annotation_to_info(variant, sample_name):
    """
    :type variant: GelVariant
    """
    """
    FILTER Fields
    """
    # Germline variant allele frequency check
    try:
        if (
            float(
                variant.annotation.additionalAttributes["GEL.GL.6628"].attribute["AF"]
            )
            >= 0.01
        ):
            if "PASS" in variant.filter:
                variant.filter.remove(u"PASS")
            variant.filter.append(u"CommonGermlineVariant")
    except (KeyError, TypeError):
        pass

    # AF_FFPE, AF_FFpcrfree, AF_FFnano somatic variant check
    try:
        for key, value in variant.annotation.additionalAttributes[
            "somatic_agg_vcf"
        ].attribute.items():
            if float(value) >= 0.05:
                if "PASS" in variant.filter:
                    variant.filter.remove(u"PASS")
                variant.filter.append(u"RecurrentSomaticVariant")
                break
    except (KeyError, TypeError):
        pass

    # High levels of sequencing noise
    try:
        fdp50 = float(variant.get_sample_value(sample_name, "FDP50")[0])
        dp50 = float(variant.get_sample_value(sample_name, "DP50")[0])
        if dp50 != 0 and fdp50 / dp50 >= 0.1:
            if "PASS" in variant.filter:
                variant.filter.remove(u"PASS")
            variant.filter.append(u"BCNoise10Indel")
    except (KeyError, StandardError):
        pass

    # Variant overlaps simple repeat
    try:
        for repeat in variant.annotation.repeat:
            if repeat.source == "trf":
                if "PASS" in variant.filter:
                    variant.filter.remove(u"PASS")
                variant.filter.append(u"SimpleRepeat")
                break
    except (KeyError, TypeError):
        pass

    # variant.annotation.populationFrequencies
    if variant.annotation.populationFrequencies:
        for pop in variant.annotation.populationFrequencies:
            if pop.study == "GNOMAD_GENOMES" and pop.population == "ALL":
                if float(pop.altAlleleFreq) >= 0.01:
                    if "PASS" in variant.filter:
                        variant.filter.remove(u"PASS")
                    variant.filter.append(u"CommonGnomADVariant")
                    break

    # Systematic Sequencing/Mapping error
    if 'SomaticFisherPhred' in variant.info:
        if float(variant.info['SomaticFisherPhred'][0]) < 50:
            if "PASS" in variant.filter:
                variant.filter.remove(u"PASS")
            variant.filter.append(u"PONnoise50SNV")




#    try:
#        fisherTest = float(variant.studies[0].files[0].attributes.SomaticFisherPhred)
#        if fisherTest < 500:
#            if "PASS" in variant.filter:
#                variant.filter.remove(u"PASS")
#            variant.filter.append(u"PONnoise50SNV")
#    except (KeyError, TypeError):
#        pass
         


    """
    INFO fields
    """
    # GNOMAD all population genome frequencies
    if variant.annotation.populationFrequencies:
        variant.info["AF_GNOMAD"] = [
            pop.altAlleleFreq
            for pop in variant.annotation.populationFrequencies
            if pop.study == "GNOMAD_GENOMES" and pop.population == "ALL"
        ]

    # Reference homopolymer intersection
    try:
        if int(variant.info["IHP"][0]) >= 8:
            variant.info["HomopolimerIndel"] = "HomopolimerIndel"
    except (KeyError, TypeError):
        pass
    # Add Segmental duplication info -> genomicSuperDup
    try:
        for repeat in variant.annotation.repeat:
            if repeat.source == "genomicSuperDup":
                variant.info["SegmentalDuplication"] = "SegmentalDuplication"
                break
    except (KeyError, TypeError):
        pass

    # GeL germline frequencies
    try:
        if variant.annotation.additionalAttributes["GEL.GL.6628"]:
            variant.info["AF_GEL_GL"] = variant.annotation.additionalAttributes[
                "GEL.GL.6628"
            ].attribute["AF"]
            variant.info["AN_GEL_GL"] = variant.annotation.additionalAttributes[
                "GEL.GL.6628"
            ].attribute["AN"]
            variant.info["AC_GEL_GL"] = variant.annotation.additionalAttributes[
                "GEL.GL.6628"
            ].attribute["AC"]
    except (KeyError, TypeError):
        pass
    # GeL somatic allele frequencies
    try:
        if variant.annotation.additionalAttributes["somatic_agg_vcf"].attribute[
            "AF_FFPE"
        ]:
            variant.info["AF_GEL_SOM_FFPE"] = variant.annotation.additionalAttributes[
                "somatic_agg_vcf"
            ].attribute["AF_FFPE"]
    except (KeyError, TypeError):
        pass
    try:
        if variant.annotation.additionalAttributes["somatic_agg_vcf"].attribute[
            "AF_FFpcrfree"
        ]:
            variant.info[
                "AF_GEL_SOM_FFpcrfree"
            ] = variant.annotation.additionalAttributes["somatic_agg_vcf"].attribute[
                "AF_FFpcrfree"
            ]
    except (KeyError, TypeError):
        pass
    try:
        if variant.annotation.additionalAttributes["somatic_agg_vcf"].attribute[
            "AF_FFnano"
        ]:
            variant.info["AF_GEL_SOM_FFnano"] = variant.annotation.additionalAttributes[
                "somatic_agg_vcf"
            ].attribute["AF_FFnano"]
    except (KeyError, TypeError):
        pass

    # Calculate VAF, extract fraction
    variant.info["VAF"] = variant.calculate_vaf(sample_name)["fraction"]

    # Format:GeneName|TranscriptID|CDSchange|ProteinChange|ConsequenceType
    consequences = []
    for cqst_type in variant.annotation.consequenceTypes:
        single_consequence = []
        cds_change_string = ""
        gene_name = cqst_type.geneName
        transcript_id = cqst_type.ensemblTranscriptId
        protein_change = ""
        sift_description = ""
        sift_score = ""
        polyphen_description = ""
        polyphen_score = ""
        if cqst_type.proteinVariantAnnotation:
            if (
                cqst_type.proteinVariantAnnotation.reference is not None
                and cqst_type.proteinVariantAnnotation.alternate is not None
            ):
                ref = cqst_type.proteinVariantAnnotation.reference
                pos = str(cqst_type.proteinVariantAnnotation.position)
                alt = [
                    x if x else "N"
                    for x in [cqst_type.proteinVariantAnnotation.alternate]
                ][0]
                protein_change = "p.(" + ref + pos + alt + ")"
            if (cqst_type.proteinVariantAnnotation.substitutionScores is not None): 
                for subScore in cqst_type.proteinVariantAnnotation.substitutionScores:
                    if subScore.source == 'sift':
                        if (subScore.description is not None):
                            sift_description = subScore.description
                            sift_description = sub(r"\s+", '_', sift_description)
                        if (subScore.score is not None):
                            sift_score = str(subScore.score)
                    elif subScore.source == 'polyphen':
                        if (subScore.description is not None):
                            polyphen_description = subScore.description
                            polyphen_description = sub(r"\s+", '_', polyphen_description)
                        if (subScore.score is not None):
                            polyphen_score = str(subScore.score)                                        
        consequence = cqst_type.sequenceOntologyTerms[0].name
        for cnqs_hgvs in variant.annotation.hgvs:
            if transcript_id is not None and transcript_id in cnqs_hgvs:
                cds_change_string = cnqs_hgvs.split(":")[1]
        single_consequence = "|".join(
            [
                x if x else ""
                for x in [
                    gene_name,
                    transcript_id,
                    cds_change_string,
                    protein_change,
                    consequence,
                    sift_description,
                    sift_score,
                    polyphen_description,
                    polyphen_score,
                ]
            ]
        )
        consequences.append(single_consequence)
    variant.info["CT"] = ",".join(consequences)
    return variant


def main(args):
    vcf_file = args.vcf_input
    annotation_file = args.json_annotation
    input_name = os.path.basename(vcf_file).split(".")[0]
    outfile_prefix = args.out_directory
    outfile = os.path.join(outfile_prefix, input_name + ".vcf.gz")
    sys.stdout.write("Starting work on: " + input_name + "\n")
    sys.stdout.write("Output directory will be: " + outfile_prefix + "\n")
    sys.stdout.write("Final path will be: " + outfile + "\n")
    # file handler for VCF file
    vcf_file_h = GVCF(vcf_file)
    # file handler for annotation file
    annotation_file_h = JsonGelVariantReader(
        annotation_file, sample_list=vcf_file_h.sample
    )
    annotation_file_h.open()

    vcf_file_h.meta_info.append(
        '##FILTER=<ID=CommonGermlineVariant,Description="Genomics England filter: variants with a population germline allele frequency above 1% in a Genomics England cohort">'
    )
    vcf_file_h.meta_info.append(
        '##FILTER=<ID=RecurrentSomaticVariant,Description="Genomics England filter: recurrent somatic variants with frequency above 5% in a Genomics England cohort">'
    )
    vcf_file_h.meta_info.append(
        '##FILTER=<ID=BCNoise10Indel,Description="Genomics England filter: average fraction of filtered basecalls within 50 bases of the indel exceeds 0.1, FDP50/DP50 > 0.1">'
    )
    vcf_file_h.meta_info.append(
        '##FILTER=<ID=SimpleRepeat,Description="Genomics England filter: variants overlapping simple repeats as defined by Tandem Repeats Finder: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz">'
    )
    vcf_file_h.meta_info.append(
        '##FILTER=<ID=CommonGnomADVariant,Description="Genomics England filter: variants with a population germline allele frequency above 1% in gnomAD dataset">'
    )
    vcf_file_h.meta_info.append(
        '##FILTER=<ID=PONnoise50SNV,Description="Genomics England filter: SomaticFisherPhred below 50, indicating somatic SNV is systematic mapping/sequencing error (applies only to SNVs that pass Strelka filters)">'
    )

    vcf_file_h.meta_info.append(
        '##FORMAT=<ID=GT, Number=.,Type=String,Description="Genotype">'
    )

    vcf_file_h.meta_info.append(
        '##INFO=<ID=CT,Number=.,Type=String,Description="Consequence type as predicted by CellBase. Format:GeneName|TranscriptID|CDSchange|ProteinChange|ConsequenceType|siftDescription|siftScore|polyphenDescription|polyphenScore">'
    )
    vcf_file_h.meta_info.append(
        '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">'
    )

    vcf_file_h.meta_info.append(
        '##INFO=<ID=AF_GNOMAD,Number=A,Type=Float,Description="Allele frequency from all populations of gnomAD genome data set">'
    )
    vcf_file_h.meta_info.append(
        '##INFO=<ID=AF_GEL_GL,Number=A,Type=Float,Description="Allele frequency from the Genomics England germline cohort">'
    )
    vcf_file_h.meta_info.append(
        '##INFO=<ID=AN_GEL_GL,Number=1,Type=Integer,Description="Total number of alleles in called genotypes from Genomics England germline cohort">'
    )
    vcf_file_h.meta_info.append(
        '##INFO=<ID=AC_GEL_GL,Number=A,Type=Integer,Description="Allele count in genotypes from Genomics England germline cohort">'
    )
    vcf_file_h.meta_info.append(
        '##INFO=<ID=AF_GEL_SOM_FFpcrfree,Number=A,Type=Float,Description="Alternate Allele Frequency in the Genomics England FFpcrfree cohort">'
    )
    vcf_file_h.meta_info.append(
        '##INFO=<ID=AF_GEL_SOM_FFnano,Number=A,Type=Float,Description="Alternate Allele Frequency in the Genomics England FFnano cohort">'
    )
    vcf_file_h.meta_info.append(
        '##INFO=<ID=AF_GEL_SOM_FFPE,Number=A,Type=Float,Description="Alternate Allele Frequency in the Genomics England FFPE cohort">'
    )
    vcf_file_h.meta_info.append(
        '##INFO=<ID=HomopolimerIndel,Number=0,Type=Flag,Description="Indels intersecting with reference homopolymers of at least 8 nucleotides">'
    )
    vcf_file_h.meta_info.append(
        '##INFO=<ID=SegmentalDuplication,Number=0,Type=Flag,Description="Variants intersecting with Segmental Duplications: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz">'
    )

    meta_info = ("\n").join(vcf_file_h.meta_info)
    header = ("\t").join(vcf_file_h.header)

    outdata = cStringIO.StringIO()
    outdata.write(meta_info + "\n")
    outdata.write(header + "\n")
    chrs = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "X", "Y", "Un", "EBV")

    for variant in annotation_file_h.read():
        if str(variant).startswith(chrs):
            variant.referenceName = "chr" + variant.referenceName
        elif str(variant).startswith("MT"):
            variant.referenceName = sub("MT", "chrM", variant.referenceName)
        elif str(variant).startswith("M"):
            variant.referenceName = "chr" + variant.referenceName
        sample_name = variant.sample_names[0]
        annotated_variant = add_annotation_to_info(variant, sample_name)
        outdata.write(str(annotated_variant) + "\n")

    p = subprocess.Popen(
        ["bcftools", "sort", "-Oz", "-o", outfile], stdin=subprocess.PIPE
    )
    p.communicate(input=outdata.getvalue())
    # with open("testing.vcf", "w") as oo:
    #     oo.write(outdata.getvalue())
    outdata.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__
    )
    parser.add_argument("-v", "--vcf_input", type=str, help="Input VCF to annotate.")
    parser.add_argument("-j", "--json_annotation", help="Input JSON to use for the annotation.")
    parser.add_argument("-o", "--out_directory", type=str, help="Output directory for annotated vcf")
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(args)
