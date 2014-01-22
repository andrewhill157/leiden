import argparse
from LeidenDatabase import *
from ExtractLeidenDataFunctions import *

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

# Process any available genes
if genes is not None:
    for gene in genes:
        print("---> " + gene + ": IN PROGRESS...")

        # Extract data and submit variants for remapping
        result = process_gene(database, gene)
        id_numbers.append(result.id_number)
        table_data.append(result.table_data)
        headers.append(result.header_data)

        if len(result.table_data) == 0:
            print(get_errors(gene))


    print('---> Retrieving Remapping Results...')
    remapping_results = get_remapping_results(id_numbers)

    print('---> Saving Output Files...')

    first_header = None

    for i in range(0, len(genes)):

            # Remember first gene with valid header row
            if first_header is None:
                first_header = i

            # Ensure consistent column order between genes
            try:
                table_data[i] = match_column_order(headers[first_header], headers[i], table_data[i])
            except ValueError:
                print('    ---> ' + genes[i] + ': Header entries do not contain expected columns.')

            # Output data to file
            output_data = format_output_text(headers[first_header], table_data[i], remapping_results[i])
            write_output_file(genes[i], output_data)

    print('---> Job Complete.')