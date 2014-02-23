from lovd.database import extract_data_functions
from lovd.database import utilities
from lovd.remapping.remapping import VariantRemapper
import os
import codecs
import argparse
import hgvs

#TODO add a description
parser = argparse.ArgumentParser(description='TODO document')

group = parser.add_mutually_exclusive_group()
parser.add_argument("file_names", help="File containing full paths to the Extracted data files to be processed.")
parser.add_argument("dna_change_column_name", help="Label for column with HGVS coordinates.")
parser.add_argument("protein_change_column_name", help="Label for column with protein change description.")
args = parser.parse_args()

with open(args.file_names, 'r') as f:
    files_to_process = f.read().splitlines()

remapper = VariantRemapper()
remapping_errors = []

for file in files_to_process:
    try:
        # Name for output file
        vcf_file_name = os.path.splitext(file)[0] + '.vcf'

        with codecs.open(file, 'r', 'utf-8') as input_file:
            file_text = [x.split('\t') for x in input_file.read().splitlines()]

            header = file_text[0]
            table_data = file_text[1:]

            hgvs_notation_index = utilities.find_string_index(header, args.dna_change_column_name)
            hgvs_notation = []

            protein_change_index = utilities.find_string_index(header, args.protein_change_column_name)
            protein_change = []

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

            vcf_file_text = extract_data_functions.format_vcf_text(vcf_notation_variants, info_column_tags)
            utilities.write_table_to_file(vcf_file_name, vcf_file_text)

    except Exception as e:
        print 'Error: ' + file + ': ' + str(e)

print 'VCF files saved adjacent files in input list.'
print 'Saving remapping errors to remapping_errors.log'
utilities.write_table_to_file('remapping_errors.log', remapping_errors)
