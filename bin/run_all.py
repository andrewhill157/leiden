"""
Andrew Hill
MacArthur Lab - 2014

Driver script for extracting and validating data from the leiden open variation databases.
"""

import subprocess
import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Driver script for extracting an validating data from the Leiden Open '
                                             'Variation Databases.')

group = parser.add_argument('-u', '--url', help='Input VCF file to be annotated.')
group = parser.add_argument('--no_download', default=False, action='store_true', help='Set to use extracted data files '
                                                                                     'rather than downloading new data.')
group = parser.add_argument('-o', '--output_directory', default='.', help='Output directory for files. Search directory '
                                                                          'for extracted data files if --use_files is '
                                                                          'specified.')
group = parser.add_argument('-f', '--force_overwrite', default=False, action='store_true', help='Overwrite existing annotation files.')
args = parser.parse_args()

output_directory = args.output_directory

if not args.no_download:
    # Download table data all genes at URL
    pipe = subprocess.Popen(['python', 'extract_data.py', '-a',
                             '-u', args.url,
                             '-o', output_directory], stdout=subprocess.PIPE)

    extract_data_log = pipe.communicate()[0]

# Produce VCF files
extracted_data_files = glob.glob(os.path.join(output_directory, '*.txt'))

files_to_annotate = []
for file in extracted_data_files:

    vcf_file_path = os.path.splitext(file)[0] + '.vcf'

    if os.path.exists(vcf_file_path) and not args.force_overwrite:
        print 'Skipping', file, ': annotated file already in output directory. Use -f to overwrite.'
    else:
        files_to_annotate.append(file)

# Write temp list of files to annotate
temp_file_name = 'files_to_annotate.temp'
with open(temp_file_name, 'w') as f:
    f.write('\n'.join(files_to_annotate))

pipe = subprocess.Popen(['python', 'generate_annotated_vcf.py', '-f', temp_file_name], stdin=subprocess.PIPE)

annotation_log = pipe.communicate()[0]

annotated_files_list = [os.path.splitext(file)[0] + '.vcf' for files in files_to_annotate]

# Validate VCF Files
with open('annotated_files.temp', 'w') as f:
    f.write(' '.join(annotated_files_list))

pipe = subprocess.Popen(['python', 'validate_annotated_vcfs.py',
                         '-f', 'annotated_files.temp',
                         '-o', os.path.join(args.output_directory, 'lovd_validated_variants.vcf'),
                         '-d', os.path.join(args.output_directory, 'lovd_discordant_variants.vcf'),
                         ], stdout=subprocess.PIPE)

validation_log = pipe.communicate()[0]