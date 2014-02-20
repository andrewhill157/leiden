import argparse
from LeidenDatabase import *
from ExtractLeidenDataFunctions import *

"""
COMMAND LINE INTERFACE
"""
parser = argparse.ArgumentParser(description="Given URL to the base URL of any lovd 2 or 3 database installation, \
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

# Get database object and print the lovd version number
print("---> DETECTING lovd VERSION: IN PROGRESS...")
database = make_leiden_database(args.leidenURL)
version_number = database.get_version_number()
print("---> DETECTING lovd VERSION: COMPLETE")
print("    ---> VERSION " + str(version_number) + " DETECTED")

id_numbers = []
table_data = []
headers = []
genes = None

# User has specified the available genes option, print a list of all available genes.
if args.availableGenes:
    print("---> CHECKING  AVAILABLE GENES...")
    print("\n".join(database.get_available_genes()))

# User has specified the all option, so extract data from all genes available on the Leiden Database
elif args.all:
    print("---> CHECKING AVAILABLE GENES...")
    genes = database.get_available_genes()

else:
    if len(args.geneID) > 0:
        genes = args.geneID
    else:
        print('Must specify at least one geneID.')


if genes is not None:
    # Process any available genes
    for gene in genes:
        print("---> " + gene + ": IN PROGRESS...")
        print("    ---> Downloading Data...")
        database.set_gene_id(gene)

        # Extract data and submit variants for remapping
        result = extract_data_and_submit_remap(database, gene)
        id_numbers.append(result.remapping_batch_id_number)
        table_data.append(result.table_entries)
        headers.append(result.column_labels)

        if len(table_data) == 0:
            print('    ---> ERROR: No Entries Found. Please Verify.')
            print('    ---> ERROR: No Entries Found. Please Verify.')


    print('---> Retrieving Remapping Results...')
    remapping_results = get_remapping_results(id_numbers)

    print('---> Saving Output Files...')

    # Output results to file for each gene
    for i in range(0, len(genes)):
        if len(remapping_results[i]) > 0:
            # Output data to file
            output_data = format_output_text(headers[i], table_data[i], remapping_results[i])

            file_name = genes[i] + '.txt'
            write_output_file(file_name, output_data)

            if args.vcf:
                # Output VCF file
                vcf_data = format_vcf_text(headers[i], table_data[i], remapping_results[i])

                file_name = genes[i] + '.vcf'
                write_output_file(file_name, vcf_data)

    print('---> Job Complete.')