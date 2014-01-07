import argparse
from LeidenDatabase import *
import traceback
import time


def save_gene_data(leiden_database, gene_id):
    """
    Given a gene_id and a valid Leiden Database URL, saves two files: <gene_id>.txt and <gene_id>_MutalyzerInput.txt.
    from the specified Leiden Database (U(http://www.dmd.nl/nmdb2/home.php?action=switch_db))

    1. <gene_id>.txt contains the extracted table data containing variants specific to the specified gene_id in the
    Leiden Database. Each variant is on its own line and columns are separated by commas. Header labels are included as
    the first line of the file.

    2. <gene_id>_MutalyzerInput.txt contains only the DNA Change column of <gene_id>.txt (one variant per line). This
    file can be directly input to the mutalyzer batch position converter tool by LOVD
    (U(https://mutalyzer.nl/batchPositionConverter))

    @param leiden_database: a LeidenDatabase object constructed using the base URL of the particular Leiden database \
    to be used.
    @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for actin,
    as specified on the Leiden Database homepage (linked above).
    """
    # Constants for file delimiters
    file_extension = '.txt'
    row_delimiter = '\n'
    column_delimiter = '\t'

    filename = "".join([gene_id, file_extension])
    mutalyzer_input_file = "".join([gene_id, "_MutalyzerInput", file_extension])

    # write table data to file in Unicode encoding (some characters are not ASCII encodable)
    with open(filename, 'w', encoding='utf-8') as f, open(mutalyzer_input_file, 'w') as mutalyzer:

        file_lines = []
        headers = leiden_database.get_table_headers(gene_id)
        hgvs_mutation_column = TextProcessing.find_string_index(headers, u'DNA\xa0change')
        headers.insert(hgvs_mutation_column + 1, u'Genomic Variant')
        file_lines.append(column_delimiter.join(headers))

        print('    ---> Extracting Data...')

        # Get table data for variants on given gene
        entries = leiden_database.get_table_data(gene_id)
        variants = [x[hgvs_mutation_column] for x in entries]

        print('    ---> Remapping Variants...')
        # Remap all variants to genomic coordinates
        remapper = VariantRemapper()

        id_number = remapper.submit_variant_batch(variants)
        

        while remapper.entries_remaining_in_batch(id_number) > 0:
            time.sleep(0.5)

        result = remapper.get_batch_results(id_number)

        f.write(row_delimiter.join(file_lines))


def get_errors(gene_id, commandline_args):
    """
    Given the list of command line arguments passed to the script, return error messages with optional stack trace.

    @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for actin,
    as specified on the Leiden Database homepage (linked above).
    @param commandline_args: a parser.parse_args() object from the argparse library. Assumed to contain an \
    argument called debug to indicate verbosity of error messages. args.debug is true, full stack traces are \
    printed for errors. A simple error message is printed otherwise. See return for details.
    @rtype: string
    @return: Error message in the format "<gene_id>: ERROR - NOT PROCESSED. Use --debug option for more details." \
    if the debug option is not specified and in the format \
    "<gene_id>: ERROR - NOT PROCESSED. STACK TRACE: \n <stack trace>" if the --debug option is specified.
    """

    # Check whether the debug option has been specified in command line arguments
    if commandline_args.debug:
        tb = traceback.format_exc()
        return "".join(["---> ", gene_id, ": ERROR - NOT PROCESSED. STACK TRACE: \n", tb])
    else:
        return "".join(["---> ", gene_id, ": ERROR - NOT PROCESSED. Use --debug option for more details."])

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

# User has specified the available genes option, print a list of all available genes.
if args.availableGenes:
    print("---> CHECKING  AVAILABLE GENES...")
    print("\n".join(database.get_available_genes()))

# User has specified the all option, so extract data from all genes available on the Leiden Database
elif args.all:
    print("---> CHECKING AVAILABLE GENES...")
    for gene in database.get_available_genes():
        try:
            print("---> " + gene + ": IN PROGRESS...")
            save_gene_data(database, gene)
            print("---> " + gene + ": COMPLETE")
        except:
            print(database.get_errors(args))

# The user has not specified all option, process their list of arguments
else:
    # No arguments passed
    if len(args.geneID) == 0:
        print("---> NO GENES PROCESSED: Must pass at least one file or use the --all option")

    # Process each gene ID the user has passed
    else:
        for gene in args.geneID:
            try:
                print("---> " + gene + ": IN PROGRESS...")
                save_gene_data(database, gene)
                print("---> " + gene + ": COMPLETE")
            except:
                print(get_errors(gene, args))
