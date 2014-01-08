import argparse
from LeidenDatabase import *
import traceback
import time


# TODO update documentation
def get_remapping_results(gene_ids, id_numbers):
    """

    @param gene_ids:
    @type gene_ids:
    @param id_numbers:
    @type id_numbers:
    @return:
    @rtype:
    """

    remapper = VariantRemapper()
    results = []

    for i in range(0, len(gene_ids)):
        print('    ---> ' + gene_ids[i])

        if id_numbers[i] > 0:
            while remapper.entries_remaining_in_batch(id_numbers[i]) > 0:
                time.sleep(0.5)

            results.append(remapper.get_batch_results(id_numbers[i]))
            print('        ---> Complete')

        else:
            results.append([''])
            print('        ---> ERROR')

    return results

# TODO update documentation
def format_output_text(header, table_data, remapping):
    """

    @param header:
    @type header:
    @param table_data:
    @type table_data:
    @param remapping:
    @type remapping:
    @return:
    @rtype:
    """
    if len(remapping) > 0:
        hgvs_mutation_column = TextProcessing.find_string_index(header, u'DNA\xa0change')
        header.insert(hgvs_mutation_column + 1, u'Genomic Variant')

        for i in range(0, len(remapping)):
            if remapping[i] is '':
                remapping[i] = 'REMAPPING_ERROR'

            table_data[i].insert(hgvs_mutation_column + 1, remapping[i])

    table_data.insert(0, header)
    return table_data


# TODO update documentation
def write_output_file(gene_id, output_data):
    """

    @param gene_id:
    @type gene_id:
    @param output_data:
    @type output_data:
    @return:
    @rtype:
    """

    # Constants for file delimiters
    file_extension = '.txt'
    row_delimiter = '\n'
    column_delimiter = '\t'

    filename = "".join([gene_id, file_extension])

    # write table data to file in Unicode encoding (some characters are not ASCII encodable)
    with open(filename, 'w', encoding='utf-8') as f:
        file_lines = []

        for row in output_data:
            file_lines.append(column_delimiter.join(row))

        f.write(row_delimiter.join(file_lines))


# TODO update documentation
def process_gene(leiden_database, gene_id):
    """

    @param leiden_database:
    @type leiden_database:
    @param gene_id:
    @type gene_id:
    @return:
    @rtype:
    """
    try:
        # Get header data
        print('    ---> Downloading Variant Data...')
        headers = leiden_database.get_table_headers(gene_id)
        hgvs_mutation_column = TextProcessing.find_string_index(headers, u'DNA\xa0change')

        print('    ---> Processing Variant Data...')
        # Get table data for variants on given gene
        entries = leiden_database.get_table_data(gene_id)
        variants = [x[hgvs_mutation_column] for x in entries]

        print('    ---> Submitting Remapping...')
        # Remap all variants to genomic coordinates
        remapper = VariantRemapper()
        id_number = remapper.submit_variant_batch(variants)

        if id_number > 0:
            print('    ---> Remapping Submitted.')
        else:
            print('    ---> ERROR: Batch Remapping Failed.')
    except:
        id_number = -1
        entries = []
        headers = []

    return {'id_number':id_number, 'table_data':entries, 'header_data':headers}


# TODO update documentation
def get_errors(gene_id):
    """

    @param gene_id:
    @return:
    @rtype:
    """

    return ''.join(['    ---> ' + 'ERROR: No Entries Found. Please Verify.'])


# TODO update documentation
"""
COMMAND LINE INTERFACE
"""
parser = argparse.ArgumentParser(description="Given a geneID, saves two files: <geneID>.txt and \
<geneID>_MutalyzerInput.txt. from the Leiden Database (http://www.dmd.nl/nmdb2/home.php?action=switch_db). \
1. <geneID>.txt contains the extracted table data containing variants specific to the specified geneID in the Leiden \
Database. Each variant is on its own line and columns are separated by commas. Header labels are included as the first \
line of the file.\
2. <geneID>_MutalyzerInput.txt contains only the DNA Change column of <geneID>.txt (one variant per line). This file \
can be directly input to the mutalyzer batch position converter tool by LOVD \
(https://mutalyzer.nl/batchPositionConverter)")
group = parser.add_mutually_exclusive_group()
group.add_argument("-d", "--debug", action="store_true",
                   help="When errors are encountered, a full stack traceback is printed.")

group2 = parser.add_mutually_exclusive_group()
group2.add_argument("-g", "--availableGenes", action="store_true", help="A list of all available genes is printed.")
group2.add_argument("-a", "--all", action="store_true",
                    help="Extract data for all available genes in the Leiden Database.")

parser.add_argument("leidenURL", help="base URL of the particular Leiden database to be used. For example, the Leiden \
 muscular dystrophy pages homepage is http://www.dmd.nl/nmdb2/. This must be a valid URL to the homepage, \
 any .php page at the end of the URL will be ignored. URLs may not contain the & character in some command line shells")
parser.add_argument("geneID", help="Gene ID or multiple geneIDs to retrieve from the Leiden Database.", nargs="*")

args = parser.parse_args()

# Get database object and print the LOVD version number
print("---> DETECTING LOVD VERSION: IN PROGRESS...")
database = get_leiden_database(args.leidenURL)
version_number = database.get_version_number()
print("---> DETECTING LOVD VERSION: COMPLETE")
print("    ---> VERSION " + str(version_number) + " DETECTED")

id_numbers = []
table_data = []
headers = []


# User has specified the available genes option, print a list of all available genes.
if args.availableGenes:
    print("---> CHECKING  AVAILABLE GENES...")
    print("\n".join(database.get_available_genes()))

# User has specified the all option, so extract data from all genes available on the Leiden Database
elif args.all:
    print("---> CHECKING AVAILABLE GENES...")
    available_genes = database.get_available_genes()
    for gene in available_genes:
        print("---> " + gene + ": IN PROGRESS...")

        # Extract data and submit variants for remapping
        result = process_gene(database, gene)
        id_numbers.append(result['id_number'])
        table_data.append(result['table_data'])
        headers.append(result['header_data'])

        if len(result['table_data']) == 0:
            print(get_errors(gene))


    print('---> Retrieving Remapping Results...')
    remapping_results = get_remapping_results(available_genes, id_numbers)

    print('---> Saving Output Files...')
    for i in range(0, len(available_genes)):
        if len(table_data[i]) > 0:
            output_data = format_output_text(headers[i], table_data[i], remapping_results[i])
            write_output_file(available_genes[i], output_data)

    print('---> Job Complete.')

# The user has not specified all option, process their list of arguments
else:
    # No arguments passed
    if len(args.geneID) == 0:
        print("---> NO GENES PROCESSED: Must pass at least one geneID or use the --all option")

    # Process each gene ID the user has passed
    else:
        for gene in args.geneID:
            print("---> " + gene + ": IN PROGRESS...")

            # Extract data and submit variants for remapping
            result = process_gene(database, gene)
            id_numbers.append(result['id_number'])
            table_data.append(result['table_data'])
            headers.append(result['header_data'])

            if len(result['table_data']) == 0:
                print(get_errors(gene))

        print('---> Retrieving Remapping Results...')
        remapping_results = get_remapping_results(args.geneID, id_numbers)

        print('---> Saving Output Files...')

        for i in range(0, len(args.geneID)):
            if len(table_data[i]) > 0:
                output_data = format_output_text(headers[i], table_data[i], remapping_results[i])
                write_output_file(args.geneID[i], output_data)

        print('---> Job Complete.')