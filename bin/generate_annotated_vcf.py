# -*- coding: utf-8 -*-
#!/usr/bin/env python

import argparse
import os
import pandas as pd
from leiden import annotate_vcf, vcf
from leiden.remapping import VariantRemapper

COLUMN_DELIMITER = '\t'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Produce an annotated VCF from raw LOVD output files. Requires that a copy'
                                                 'of variant_effect_predictor is on PATH with cache 27 and 28 installed and'
                                                 'the Downstream plugin installed.')

    group = parser.add_argument('-i', '--input_files',  help='List of input raw LOVD output files to be annotated.', nargs='*')

    args = parser.parse_args()

    rm = VariantRemapper()

    for file in args.input_files:
        base_file_name = os.path.splitext(file)[0]
        annotation_input_file = base_file_name + '_HGVS.temp'
        output_file = base_file_name + '.vcf'

        # Clean LOVD data for VCF
        lovd_file = pd.read_csv(file, sep=COLUMN_DELIMITER)
        lovd_file = vcf.remove_malformed_fields(lovd_file)

        # Output VCF variants for annotation
        column_list = ['dna_change', 'protein_change', 'var_pub_as', 'rna_change', 'db_id', 'variant_remarks', 'reference', 'frequency']


        vcf_format = vcf.convert_to_vcf_format(lovd_file[column_list], rm, 'dna_change', 'LOVD')
        vcf_column_order = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

        with open(annotation_input_file, 'w') as f:
            vcf_header = ['##fileformat=VCFv4.0',
                          vcf.get_vcf_info_header(lovd_file[column_list], 'LOVD', 'Data from LOVD'),
                          '#' + '\t'.join(vcf_column_order)
            ]
            f.write('\n'.join(vcf_header) + '\n')

            vcf_format[vcf_column_order].to_csv(f, sep=COLUMN_DELIMITER, header=False, index=False)

        # Annotate with VEP
        annotate_vcf.annotate_vep(annotation_input_file, output_file)

        # Combine VEP VCF output with original LOVD data
        with open(output_file, 'r') as f:
            header_lines = vcf.get_vcf_header_lines(f)

        cumulative_vcf = pd.read_csv(output_file, header=len(header_lines)-1, sep=COLUMN_DELIMITER)

        # Make sure only variants with both CSQ and LOVD tags are in INFO column (VEP can't annotate some)
        cumulative_vcf = cumulative_vcf[cumulative_vcf['INFO'].str.contains('CSQ') & cumulative_vcf['INFO'].str.contains('LOVD')]

        # Write final VCF
        with open(output_file, 'w') as f:
            f.writelines('\n'.join(header_lines) + '\n')
            cumulative_vcf.to_csv(f, header=False, index=False, sep=COLUMN_DELIMITER)

        os.remove(annotation_input_file)