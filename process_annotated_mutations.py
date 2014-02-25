"""
Andrew Hill
MacArthur Lab - 2014

Outputs error logs and validation statistics on VEP and 26K annotated VCF files. Meant to validate variants converted
from HGVS to VCF format.

For help, execute: python process_annotated_mutations.py --help
"""

import argparse
import os
from lovd.io.file_io import write_table_to_file
from lovd.validation import annotation_processing
from lovd.io import file_io

parser = argparse.ArgumentParser(description='Output validation statitics and error logs when processing VEP and AC '
                                             'annotated VCF files. Outputs processing_errors.log and '
                                             'discordant_annotations.log, which contain information on variants that '
                                             'could not be validated due to syntax errors and variants that not '
                                             'correctly validate respectively.')

group = parser.add_argument_group()
group.add_argument('-f', '--file_names', required=True, help="File containing full paths to the VCF files to be processed.")
group.add_argument('-o', '--output_directory', default='.', help='Output directory for .log files.')
args = parser.parse_args()

# Make the output directory if does not already exist
output_directory = args.output_directory

if not os.path.exists(output_directory):
    os.mkdir(args.output_directory)


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
    vcf_file = file_io.read_table_from_file(file)
    variants = [x for x in vcf_file if not x[0].startswith('#')]

    for variant_row in variants:
        # Extract info from line of VCF
        chromosome_number = variant_row[0]
        coordinate = variant_row[1]
        ref = variant_row[3]
        alt = variant_row[4]
        info = variant_row[7].split(';')

        # Get link to relevant location in UCSC genome browser
        viewing_interval = 25
        if coordinate != '.':
            ucsc_link = annotation_processing.get_ucsc_location_link(chromosome_number,
                                                                     str(int(coordinate) - viewing_interval),
                                                                     str(int(coordinate) + viewing_interval))

        # Get tagged information from the info column
        severe_impact = annotation_processing.get_tagged_entry_value(info, 'SEVERE_IMPACT')

        try:
            hgvs_mutation = annotation_processing.get_tagged_entry_value(info, 'HGVS')
        except:
            hgvs_mutation = 'NOT_FOUND'

        try:
            protein_change = annotation_processing.get_tagged_entry_value(info, 'LAA_CHANGE')
        except:
            protein_change = 'NOT_FOUND'

        # Check concordance
        try:
            laa_change = ''
            aa_change = ''

            laa_change = annotation_processing.get_laa_change(info)
            aa_change = annotation_processing.get_aa_change(info)

            if annotation_processing.is_concordant_annotation(laa_change, aa_change):
                concordant_annotation_count += 1
            else:
                discordant_annotation_count += 1
                discordant_annotations.append([file,
                                               hgvs_mutation,
                                               chromosome_number,
                                               coordinate,
                                               ref,
                                               alt,
                                               ucsc_link,
                                               protein_change,
                                               annotation_processing.map_aa_codes(laa_change.before) + '/' +
                                               annotation_processing.map_aa_codes(laa_change.after),
                                               annotation_processing.map_aa_codes(aa_change.before) + '/' +
                                               annotation_processing.map_aa_codes(aa_change.after), severe_impact])
        except Exception as e:
            error_count += 1
            processing_errors.append([file,
                                      str(e),
                                      hgvs_mutation,
                                      protein_change])

        # Check allele frequency
        allele_frequency = annotation_processing.get_overall_26K_allele_frequency(info)

        if allele_frequency > 0.5:
            high_26K_frequency_count += 1
        if allele_frequency > 0:
            overlap_26K_count += 1

        # Check HGMD Overlap
        HGMD_SITE_TAG = 'HGMD_SITE'
        HGMD_MUTATION_TAG = 'HGMD_MUT'
        hgmd_site = annotation_processing.get_tagged_entry_value(info, 'HGMD_SITE')
        hgmd_mutation = annotation_processing.get_tagged_entry_value(info, 'HGMD_MUT')

        if hgmd_site != '':
            hgmd_site_count += 1
        if hgmd_mutation != '':
            hgmd_mutation_count += 1

        total_mutation_count += 1

# Write out errors and discordant mutations to file
processing_errors.insert(0, ['file',
                             'error',
                             'hgvs',
                             'protein'])

processing_errors_file_name = os.path.join(output_directory, 'processing_errors.log')
write_table_to_file(processing_errors_file_name, processing_errors)

discordant_annotations.insert(0, ['file',
                                  'hgvs',
                                  'chromosome',
                                  'coordinate',
                                  'ref',
                                  'alt',
                                  'ucsc',
                                  'protein',
                                  'aa_change_lovd',
                                  'aa_change_vep',
                                  'severe_impact'])
discordant_annotations_file_name = os.path.join(output_directory, 'discordant_annotations.log')
write_table_to_file(discordant_annotations_file_name, discordant_annotations)

# Print results
print 'Total Mutations: ' + str(total_mutation_count)
print 'Concordant Annotations: ' + str(concordant_annotation_count)
print 'Discordant Annotations: ' + str(discordant_annotation_count)
print 'Errors and Indels: ' + str(error_count)
print ''
print 'HGMD Sites: ' + str(hgmd_site_count)
print 'HGMD Mutations: ' + str(hgmd_mutation_count)
print ''
print 'High 26K Frequency: ' + str(high_26K_frequency_count)
print '26K Overlap Count: ' + str(overlap_26K_count)












    


            