"""
Andrew Hill
MacArthur Lab - 2014

Script that makes makes use of LOVD modules to generate VCF files from files containing HGVS variants and their
predicted protein changes.

For help, execute: python generate_vcf_files.py --help
"""

import os
import argparse

from lovd.database import utilities
from lovd.io import file_io
from lovd.remapping.remapping import VariantRemapper


#TODO add a description
parser = argparse.ArgumentParser(description='Generates VCF files from a tab-separated file containing a list of '
                                             'mutations in HGVS notation and their corresponding predicted protein '
                                             'change (for validation purposes). The first line of the file must contain '
                                             'labels for the columns. Must include the label used for each column as '
                                             'command line parameters. Defaults are dna and protein. '
                                             'This script can be used to regenerate VCF files from extracted data from '
                                             'the leiden databases if needed, or generate a VCF using HGVS data from '
                                             'any other source. '
                                             'Errors from remapping coordinates are listed in remapping_errors.log in '
                                             'the script directory or in output directory if specified. '
                                             'VCF files are generated adjacent to their source or in output directory '
                                             'if specified. Output file does not include file mutations that failed to remap.')

group = parser.add_argument_group()
group.add_argument('-l', '--file_list', required=True, help='File containing list of full paths to the data files to be processed.')
group.add_argument('-m', '--hgvs_mutation_column', required=True, help='Label for column with HGVS coordinates.')
group.add_argument('-p', '--protein_change_column', required=True, help='Label for column with protein change description.')
group.add_argument('-o', '--output_directory', default='', help='Output directory for saved files.')
args = parser.parse_args()

# Make the output directory if does not already exist
output_directory = args.output_directory

if output_directory != '' and not os.path.exists(output_directory):
    os.mkdir(args.output_directory)


# Read list of files to process
with open(args.file_list, 'r') as f:
    files_to_process = f.read().splitlines()

remapper = VariantRemapper()
remapping_errors = []

for file in files_to_process:
    try:
        table_data = file_io.read_table_from_file(file)
        header = table_data.pop(0)

        # Isolate data columns with HGVS mutations and protein change
        hgvs_notation_index = utilities.find_string_index(header, args.hgvs_mutation_column)
        hgvs_notation = []

        protein_change_index = utilities.find_string_index(header, args.protein_change_column)
        protein_change = []

        # Remap Variants and build list for INFO column tags
        vcf_notation_variants = []

        for row in table_data:
            try:
                vcf_notation = remapper.hgvs_to_vcf(row[hgvs_notation_index])
                vcf_notation_variants.append(vcf_notation)

                hgvs_notation.append(row[hgvs_notation_index])
                protein_change.append(row[protein_change_index])

            except Exception as e:
                remapping_errors.append([file, row[hgvs_notation_index], str(e)])

        info_column_tags = {'HGVS': ('string', 'LOVD HGVS notation describing DNA change', hgvs_notation),
                'LAA_CHANGE': ('string', 'LOVD amino acid change', protein_change)}


        if output_directory != '':
            # Save output files to specified directory
            base_file_name = os.path.basename(file)
        else:
            # Save output files adjacent to input files
            base_file_name = file

        # Construct VCF file path
        vcf_file_name = os.path.splitext(base_file_name)[0] + '.vcf'
        vcf_file_name = os.path.join(output_directory, vcf_file_name)

        vcf_file_text = file_io.format_vcf_text(vcf_notation_variants, info_column_tags)
        file_io.write_table_to_file(vcf_file_name, vcf_file_text)

    # Display any files that failed to process entirely
    except Exception as e:
        print 'Error: ' + file + ': ' + str(e)

# Remind user of file locations and save error log
error_log_file_name = os.path.join(output_directory, 'remapping_errors.log')
file_io.write_table_to_file(error_log_file_name, remapping_errors)
print 'Remapping errors to: ' + error_log_file_name
