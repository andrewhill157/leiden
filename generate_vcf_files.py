"""
Andrew Hill
MacArthur Lab - 2014

Script that makes makes use of LOVD modules to generate VCF files from files containing HGVS variants and their
predicted protein changes.

For help, execute: python generate_vcf_files.py --help
"""

import os
import argparse
from leiden.io import file_io

from lovd import utilities
from remapping.remapping import VariantRemapper


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
group.add_argument('-f', '--file_list', required=True, help='File containing list of full paths to the data files to be processed.')
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
        lovd_hgvs_notation_index = utilities.find_string_index(header, 'DNA')
        lovd_hgvs_notation = []

        lovd_db_id_index = utilities.find_string_index(header, 'DB-ID')
        lovd_db_id = []

        lovd_genetic_origin_index = utilities.find_string_index(header, 'GENET_ORI')
        lovd_genetic_origin = []

        lovd_reference_index = utilities.find_string_index(header, 'REFERENCE')
        lovd_reference = []

        lovd_template_index = utilities.find_string_index(header, 'TEMPLATE')
        lovd_template = []

        lovd_technique_index = utilities.find_string_index(header, 'TECHNIQUE')
        lovd_technique = []

        protein_change_index = utilities.find_string_index(header, 'PROTEIN')
        protein_change = []

        # Remap Variants and build list for INFO column tags
        vcf_notation_variants = []

        for row in table_data:

            lovd_hgvs_notation.append(row[lovd_hgvs_notation_index].replace(r'\s', ''))
            lovd_db_id.append(row[lovd_db_id_index].replace(r'\s', ''))
            lovd_genetic_origin.append(row[lovd_genetic_origin_index].replace(r'\s', ''))
            lovd_reference.append(row[lovd_reference_index].replace(r'\s', ''))
            lovd_template.append(row[lovd_template_index].replace(r'\s', ''))
            lovd_technique.append(row[lovd_technique_index].replace(r'\s', ''))
            protein_change.append(row[protein_change_index].replace(r'\s', ''))

            try:
                vcf_notation = remapper.hgvs_to_vcf(row[lovd_hgvs_notation_index])
                vcf_notation_variants.append(vcf_notation)
            except Exception as e:
                remapping_errors.append([file, row[lovd_hgvs_notation_index], str(e)])

        info_column_tags = {'LOVD_HGVS': ('String', 'LOVD HGVS notation describing DNA change', lovd_hgvs_notation),
                            'LOVD_DB_ID': ('String', 'LOVD database ID', lovd_db_id),
                            'LOVD_GENETIC_ORIGIN': ('String', 'LOVD genetic origin entry', lovd_genetic_origin),
                            'LOVD_REFERENCE': ('String', 'LOVD publication references', lovd_reference),
                            'LOVD_TEMPLATE': ('String', 'LOVD template entry', lovd_template),
                            'LOVD_TECHNIQUE': ('String', 'LOVD techniques entry', lovd_technique),
                            'LOVD_AA_CHANGE': ('String', 'LOVD amino acid change', protein_change)}

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
