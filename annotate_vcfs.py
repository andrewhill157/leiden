"""
Andrew Hill
MacArthur Lab - 2014

Script that makes use of the macarthur_core modules to annotate a given VCF file with 26K, HGMD, and VEP. Note that
this script can only be run on the Broad Cluster and is dependent on absolute paths to given

For help, execute: python annotate_vcfs.py --help
"""

import argparse

# Command line interface definition
from leiden.broad_cluster import annotate_vcf

parser = argparse.ArgumentParser(description='Annotate a given file using 26K, HGMD, and VEP. Note that this script '
                                             'requires the macarthur_core package and must be run on the Broad Cluster, '
                                             'as it depends on absolute paths to tools housed there.')

group = parser.add_argument('-i', '--input_file', help='Input VCF file to be annotated.')
group = parser.add_argument('-o', '--output_file', help='Output file for annotated VCF.')

args = parser.parse_args()

# Annotate with VEP and 26K data
annotate_vcf.annotate_vep(args.input_file, args.output_file)
annotate_vcf.annotate_26k(args.output_file, args.output_file)
annotate_vcf.annotate_hgmd(args.output_file, args.output_file)
annotate_vcf.annotate_dbsnp(args.output_file, args.output_file)
