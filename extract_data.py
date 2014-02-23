import argparse
from lovd.database.extract_data_functions import *
from lovd.database.utilities import write_table_to_file
from lovd.database.leiden_database import make_leiden_database
from lovd.remapping.remapping import VariantRemapper

"""
COMMAND LINE INTERFACE
"""
parser = argparse.ArgumentParser(description="Given URL to the base URL of any LOVD 2 or 3 database installation, \
    such as http://www.dmd.nl/nmdb2/, extract variant entries associated with specified genes. Can specify either \
    a space-separated list of gene names or -a option to extract data from all genes at the specified URL. Variants \
    for each gene are saved in a file named according to the gene name they are associated with. Optional VCF output \
    also available via -v option")

group1 = parser.add_mutually_exclusive_group()
group1 = group1.add_argument('-v', '--vcf', action='store_true', help='Output VCF formatted files for each gene. Regular data \
                                                               files will still be output. VCF files named according \
                                                               to <geneID>_VCF.txt')

group2 = parser.add_mutually_exclusive_group()
group2.add_argument("-g", "--availableGenes", action="store_true", help="A list of all available genes is printed.")
group2.add_argument("-a", "--all", action="store_true",
                    help="Extract data for all available genes in the specified Leiden Database.")

parser.add_argument("leidenURL", help="base URL of the particular Leiden database to be used. For example, the Leiden \
 muscular dystrophy pages homepage is http://www.dmd.nl/nmdb2/. This must be a valid URL to base page of database. For \
 example, http://databases.lovd.nl/whole_genome/ is a valid LOVD3 URL, while http://databases.lovd.nl/whole_genome/genes \
 is not. A list of such acceptable URLs is maintained here: http://www.lovd.nl/2.0/index_list.php")
parser.add_argument("geneID", help="Gene ID or multiple geneIDs to retrieve from the Leiden Database.", nargs="*")

args = parser.parse_args()

genes = None

# User has specified the available genes option, print a list of all available genes.
if args.availableGenes:
    database = make_leiden_database(args.leidenURL)

    print("\n".join(database.get_available_genes()))

else:
    # Get database object and print the lovd version number
    print("---> DETECTING LOVD VERSION: IN PROGRESS...")
    database = make_leiden_database(args.leidenURL)
    version_number = database.get_version_number()
    print("---> DETECTING LOVD VERSION: COMPLETE")
    print("    ---> VERSION " + str(version_number) + " DETECTED")

    # User has specified the all option, so extract data from all genes available on the Leiden Database
    if args.all:
        print("---> CHECKING AVAILABLE GENES...")
        genes = database.get_available_genes()

    else:
        if len(args.geneID) > 0:
            genes = args.geneID
        else:
            print('Must specify at least one geneID.')


    if genes is not None:

        print('---> Loading refSeq transcripts for remapping...')
        remapper = VariantRemapper()

        remapping_errors = []
        variant_entries_found = False

        for gene in genes:
            print '---> ' + gene + ': IN PROGRESS...'
            print '    ---> Downloading data...'

            try:
                table_data, column_labels = extract_data(database, gene)
                variant_entries_found = True
            except ValueError as e:
                print '    ---> ' + str(e)
                variant_entries_found = False

            if variant_entries_found:
                hgvs_index = utilities.find_string_index(column_labels, 'dna')
                protein_change_index = utilities.find_string_index(column_labels, 'protein')

                print '    ---> Remapping variants'
                vcf_format_variants = []
                hgvs_notation = []
                protein_change = []

                for row in table_data:
                    try:
                        remapped_variant = remapper.hgvs_to_vcf(row[hgvs_index])
                        hgvs_notation.append(row[hgvs_index])

                        protein_change.append(row[protein_change_index])

                        vcf_format_variants.append(remapped_variant)

                    except Exception as e:
                        remapping_errors.append([gene, row[hgvs_index], str(e)])

                info_column_tags = {'HGVS': ('string', 'LOVD HGVS notation describing DNA change', hgvs_notation),
                                    'LAA_CHANGE': ('string', 'LOVD amino acid change', protein_change)}

                vcf_text = format_vcf_text(vcf_format_variants, info_column_tags)

                print '    ---> Saving VCF file...'
                vcf_file_name = gene + '.vcf'
                utilities.write_table_to_file(vcf_file_name, vcf_text)

                print '    ---> Saving raw data...'
                output_file_name = gene + '.txt'
                table_data.insert(0, column_labels)

                utilities.write_table_to_file(output_file_name, table_data)

        print '---> Saving remapping errors to remapping_errors.log'
        utilities.write_table_to_file('remapping_errors.log', remapping_errors)

        print('---> All genes complete..')