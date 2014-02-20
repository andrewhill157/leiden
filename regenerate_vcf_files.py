from lovd.database import extract_data_functions
from lovd.database.utilities import Utilities
from collections import namedtuple
import os
import argparse

#TODO add a description
parser = argparse.ArgumentParser(description='TODO document')

group = parser.add_mutually_exclusive_group()
parser.add_argument("file_names", help="File containing full paths to the Extracted data files to be processed.")
args = parser.parse_args()

with open(args.file_names, 'r') as f:
    files_to_process = f.read().splitlines()

for file in files_to_process:
    print(file)
    try:
        # Name for output file
        vcf_file_name = os.path.splitext(file)[0] + '.vcf'

        with open(file, 'r') as input_file:
            file_text = [x.split('\t') for x in input_file.read().splitlines()]

            header = file_text[0]
            table_data = file_text[1:]

            chromosome_number = [x[Utilities.find_string_index(header, 'Chromosome number')] for x in table_data]
            coordinate = [x[Utilities.find_string_index(header, 'coordinate')] for x in table_data]
            ref = [x[Utilities.find_string_index(header, 'ref')] for x in table_data]
            alt = [x[Utilities.find_string_index(header, 'alt')] for x in table_data]

            remapping_results = namedtuple('remapping_results', 'chromosome_number coordinate ref alt')
            remapping = remapping_results(chromosome_number, coordinate, ref, alt)

            vcf_text = extract_data_functions.format_vcf_text(header, table_data, remapping)
            extract_data_functions.write_output_file(vcf_file_name, vcf_text)
    except Exception as e:
        print('Error: ' + file)
        for lines in table_data:
            print(str(e) + file)