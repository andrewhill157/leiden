"""
Andrew Hill 
MacArthur Lab - 2014

Script that makes use of core LOVD module to extract variant data for different genes on any LOVD 2 or 3 database.
Also capable of remapping variants and saving to VCF files via --vcf option.

For help, execute: python extract_data.py --help
"""

import argparse
from collections import namedtuple
import os

from lovd.io import file_io
from lovd.database import utilities
from lovd.database.leiden_database import make_leiden_database
from lovd.remapping.remapping import VariantRemapper


# Command line interface definition
parser = argparse.ArgumentParser(description='Given URL to the base URL of any LOVD 2 or 3 database installation, '
                                             'such as http://www.dmd.nl/nmdb2/, extract variant entries associated with '
                                             'specified genes. Can specify either a space-separated list of gene names '
                                             'or -a option to extract data from all genes at the specified URL. Variants '
                                             'for each gene are saved in a file named according to the gene name they '
                                             'are associated with. Optional VCF output also available via -v option')

group1 = parser.add_mutually_exclusive_group()
group1 = group1.add_argument('-v', '--vcf', action='store_true', help='Output VCF formatted files for each gene. Regular data \
                                                               files will still be output. VCF files named according \
                                                               to <gene_list>_VCF.txt')

group2 = parser.add_mutually_exclusive_group()
group2.add_argument('-g', '--genes_available', action="store_true", help='Set to list of all available genes is printed. '
                                                                        'No data is extracted or written to file.')
group2.add_argument('-a', '--all', default=False, action='store_true',
                    help='Set to extract data for all available genes in the specified Leiden Database.')

parser.add_argument('-u', '--leiden_url', required=True, help='base URL of the particular Leiden database to be used. For example, '
                                               'the Leiden muscular dystrophy pages homepage is http://www.dmd.nl/nmdb2/. '
                                               'This must be a valid URL to base page of database. For example, '
                                               'http://databases.lovd.nl/whole_genome/ is a valid LOVD3 URL, while '
                                               'http://databases.lovd.nl/whole_genome/genes is not (it is not the base). '
                                               'A list of such acceptable URLs is maintained here: '
                                               'http://www.lovd.nl/2.0/index_list.php')

parser.add_argument('-l', '--gene_list', help='Gene ID or multiple gene_lists to retrieve from the Leiden Database.', nargs='*')

parser.add_argument('-o', '--output_directory', default='.', help='Output directory for saved files.')
args = parser.parse_args()

# Make the output directory if does not already exist
output_directory = args.output_directory

if not os.path.exists(output_directory):
    os.mkdir(args.output_directory)

genes = None

# User has specified the available genes option, print a list of all available genes.
def extract_data(leiden_database, gene_id):
    """
    Extracts variant table data for given gene in leiden_database and submits variants for remapping to genomic coordinates.

    @param leiden_database: database containing tables of variant data for specified gene_id
    @type leiden_database: LeidenDatabase
    @param gene_id: a string with the Gene ID of the gene to be extracted.
    @type gene_id: string
    @return: namedtuple containing remapping_batch_id_number, table_entries, and column_label entries.
    @rtype: namedtuple
    @raise: IOError if could not get data
    """
    try:
        leiden_database.set_gene_id(gene_id)
        column_labels = leiden_database.get_table_headers()
        table_entries = leiden_database.get_table_data()

    except Exception as e:
        raise e

    results = namedtuple('results', 'table_entries, column_labels')
    return results(table_entries, column_labels)


if args.genes_available:
    database = make_leiden_database(args.leiden_url)

    print("\n".join(database.get_available_genes()))

else:
    # Get database object and print the lovd version number
    print("---> DETECTING LOVD VERSION: IN PROGRESS...")
    
    # Use factory method to get LeidenDatabase object
    database = make_leiden_database(args.leiden_url)
    version_number = database.get_version_number()
    
    print("---> DETECTING LOVD VERSION: COMPLETE")
    print("    ---> VERSION " + str(version_number) + " DETECTED")

    # User has specified the all option, so extract data from all genes available on the Leiden Database
    if args.all:
        print("---> CHECKING AVAILABLE GENES...")
        genes = database.get_available_genes()

    else:
        if len(args.gene_list) > 0:
            genes = args.gene_list
        else:
            print('Must specify at least one gene_list.')


    if genes is not None:

        print('---> Loading refSeq transcripts for remapping...')
        remapper = VariantRemapper()

        remapping_errors = []
        variant_entries_found = False

        # Extract data and output VCF and raw data for each gene
        for gene in genes:
            print '---> ' + gene + ': IN PROGRESS...'
            print '    ---> Downloading data...'

            # Extract table data
            try:
                table_data, column_labels = extract_data(database, gene)
                variant_entries_found = True
            except ValueError as e:
                print '    ---> ' + str(e)
                variant_entries_found = False

            # Write output files only if variants entries are found
            if variant_entries_found:

                # Generate VCF files if option specified
                if args.vcf:
                    # Extract info for INFO column tags
                    hgvs_index = utilities.find_string_index(column_labels, 'dna')
                    protein_change_index = utilities.find_string_index(column_labels, 'protein')

                    print '    ---> Remapping variants'
                    vcf_format_variants = []
                    hgvs_notation = []
                    protein_change = []

                    # Remap all variants and form list for INFO column tags
                    for row in table_data:
                        try:
                            remapped_variant = remapper.hgvs_to_vcf(row[hgvs_index])
                            hgvs_notation.append(row[hgvs_index])

                            protein_change.append(row[protein_change_index])

                            vcf_format_variants.append(remapped_variant)

                        except Exception as e:
                            # Maintain a list of errors that occur during remapping
                            remapping_errors.append([gene, row[hgvs_index], str(e)])

                    info_column_tags = {'HGVS': ('string', 'LOVD HGVS notation describing DNA change', hgvs_notation),
                                        'LAA_CHANGE': ('string', 'LOVD amino acid change', protein_change)}

                    vcf_text = file_io.format_vcf_text(vcf_format_variants, info_column_tags)

                    print '    ---> Saving VCF file...'
                    vcf_file_name = os.path.join(output_directory, gene + '.vcf')
                    file_io.write_table_to_file(os.path.join(vcf_file_name), vcf_text)

                print '    ---> Saving raw data...'
                table_data.insert(0, column_labels)
                output_file_name = os.path.join(output_directory, gene + '.txt')

                file_io.write_table_to_file(output_file_name, table_data)

        print '---> Saving remapping errors to remapping_errors.log'
        error_log_file_name = os.path.join(output_directory, 'remapping_errors.log')
        file_io.write_table_to_file(error_log_file_name, remapping_errors)

        print('---> All genes complete..')