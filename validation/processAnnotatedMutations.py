import argparse
import AnnotationProcessing
import logging

ERROR_FILENAME = 'errors.log'
DISCORDANT_FILENAME = 'discordant.log'
logging.basicConfig(filename=ERROR_FILENAME,filemode='w',level=logging.DEBUG)
logging.basicConfig(filename=DISCORDANT_FILENAME,filemode='w',level=logging.WARNING)

#TODO add a description
parser = argparse.ArgumentParser(description='TODO document')

group = parser.add_mutually_exclusive_group()
parser.add_argument("file_names", help="File containing full paths to the VCF files to be processed.")
args = parser.parse_args()

# Counts of various categories of variants
concordant_annotation_count = 0
discordant_annotation_count = 0
error_count = 0

hgmd_site_count = 0
hgmd_mutation_count = 0

high_26K_frequency_count = 0
overlap_26K_count = 0

total_mutation_count = 0

with open(args.file_names, 'r') as file_list:
    files_to_process = file_list.read().splitlines()

for file in files_to_process:
    # Extract all non-header lines from VCF file
    with open(file, 'r') as vcf_file:
        variants = [x for x in vcf_file if not x.startswith('#')]

    for variant in variants:

        vcf_info_column_list = AnnotationProcessing.get_vcf_info_column(variant)

        severe_impact = AnnotationProcessing.get_severe_impact(vcf_info_column_list)

        # Get original HGVS notation
        try:
            hgvs_mutation = AnnotationProcessing.get_tagged_entry_value(vcf_info_column_list, 'HGVS')
        except:
            hgvs_mutation = 'NOT_FOUND'

        # Get original protein change
        try:
            protein_change = AnnotationProcessing.get_tagged_entry_value(vcf_info_column_list, 'LAA_CHANGE')
        except:
            protein_change = 'NOT_FOUND'

        # Check concordance
        try:
            laa_change = ''
            aa_change = ''

            laa_change = AnnotationProcessing.get_laa_change(vcf_info_column_list)
            aa_change = AnnotationProcessing.get_aa_change(vcf_info_column_list)

            if AnnotationProcessing.is_concordant_annotation(laa_change, aa_change):
                concordant_annotation_count += 1
            else:
                discordant_annotation_count += 1
                logging.warning('FILE_NAME: ' + file + '; HGVS:' + hgvs_mutation + '; PROTEIN_CHANGE: ' + protein_change +
                                ';LAA_CHANGE: ' + str(laa_change) + '; AA_CHANGE: ' + str(aa_change) + '; SEVERE_IMPACT: ' + severe_impact)
        except Exception as e:
            logging.debug('FILE_NAME: ' + file + '; ERROR_MESSAGE: ' + str(e) + '; HGVS: ' + hgvs_mutation + '; PROTEIN_CHANGE: ' + protein_change)
            error_count += 1

        # Check allele frequency
        allele_frequency = AnnotationProcessing.get_overall_26K_allele_frequency(vcf_info_column_list)

        if allele_frequency > 0.5:
            high_26K_frequency_count += 1
        if allele_frequency > 0:
            overlap_26K_count += 1

        # Check HGMD Overlap
        HGMD_SITE_TAG = 'HGMD_SITE'
        HGMD_MUTATION_TAG = 'HGMD_MUT'
        hgmd_site = AnnotationProcessing.get_tagged_entry_value(vcf_info_column_list, 'HGMD_SITE')
        hgmd_mutation = AnnotationProcessing.get_tagged_entry_value(vcf_info_column_list, 'HGMD_MUT')

        if hgmd_site != '':
            hgmd_site_count += 1
        if hgmd_mutation != '':
            hgmd_mutation_count += 1

        total_mutation_count += 1


# Print results
print('Total Mutations: ' + str(total_mutation_count))
print('Concordant Annotations: ' + str(concordant_annotation_count))
print('Discordant Annotations: ' + str(discordant_annotation_count))
print('Errors and Indels: ' + str(error_count))
print('')
print('HGMD Sites: ' + str(hgmd_site_count))
print('HGMD Mutations: ' + str(hgmd_mutation_count))
print('')
print('High 26K Frequency: ' + str(high_26K_frequency_count))
print('26K Overlap Count: ' + str(overlap_26K_count))












    


            