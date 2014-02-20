import argparse
from lovd.database import utilities
from lovd.validation import annotation_processing

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

processing_errors = []
discordant_annotations = []

with open(args.file_names, 'r') as file_list:
    files_to_process = file_list.read().splitlines()

for file in files_to_process:
    # Extract all non-header lines from VCF file
    with open(file, 'r') as vcf_file:
        variants = [x for x in vcf_file if not x.startswith('#')]

    for variant in variants:

        vcf_info_column_list = annotation_processing.get_vcf_info_column(variant)

        severe_impact = annotation_processing.get_severe_impact(vcf_info_column_list)

        # Get original HGVS notation
        try:
            hgvs_mutation = annotation_processing.get_tagged_entry_value(vcf_info_column_list, 'HGVS')
        except:
            hgvs_mutation = 'NOT_FOUND'

        # Get original protein change
        try:
            protein_change = annotation_processing.get_tagged_entry_value(vcf_info_column_list, 'LAA_CHANGE')
        except:
            protein_change = 'NOT_FOUND'

        # Check concordance
        try:
            laa_change = ''
            aa_change = ''

            laa_change = annotation_processing.get_laa_change(vcf_info_column_list)
            aa_change = annotation_processing.get_aa_change(vcf_info_column_list)

            if annotation_processing.is_concordant_annotation(laa_change, aa_change):
                concordant_annotation_count += 1
            else:
                discordant_annotation_count += 1
                discordant_annotations.append([file, hgvs_mutation, protein_change,
                                               annotation_processing.map_aa_codes(laa_change.before) + '/' +
                                               annotation_processing.map_aa_codes(laa_change.after),
                                               annotation_processing.map_aa_codes(aa_change.before) + '/' +
                                               annotation_processing.map_aa_codes(aa_change.after), severe_impact])
        except Exception as e:
            error_count += 1
            processing_errors.append([file, str(e), hgvs_mutation, protein_change])

        # Check allele frequency
        allele_frequency = annotation_processing.get_overall_26K_allele_frequency(vcf_info_column_list)

        if allele_frequency > 0.5:
            high_26K_frequency_count += 1
        if allele_frequency > 0:
            overlap_26K_count += 1

        # Check HGMD Overlap
        HGMD_SITE_TAG = 'HGMD_SITE'
        HGMD_MUTATION_TAG = 'HGMD_MUT'
        hgmd_site = annotation_processing.get_tagged_entry_value(vcf_info_column_list, 'HGMD_SITE')
        hgmd_mutation = annotation_processing.get_tagged_entry_value(vcf_info_column_list, 'HGMD_MUT')

        if hgmd_site != '':
            hgmd_site_count += 1
        if hgmd_mutation != '':
            hgmd_mutation_count += 1

        total_mutation_count += 1

# Write out errors and discordant mutations to file
processing_errors.insert(0, ['file', 'error', 'hgvs', 'protein'])
utilities.write_output_file('processing_errors.log', processing_errors)

discordant_annotations.insert(0, ['file', 'hgvs', 'protein', 'aa_change_lovd', 'aa_change_vep', 'severe_impact'])
utilities.write_output_file('discordant_annotations.log', discordant_annotations)

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












    


            