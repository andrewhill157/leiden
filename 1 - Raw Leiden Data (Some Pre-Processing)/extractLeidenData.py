import argparse
from LeidenDatabase import *


parser = argparse.ArgumentParser(description="Given a geneID, saves two files: <geneID>.txt and \
<geneID>_MutalizerInput.txt. from the Leiden Database (http://www.dmd.nl/nmdb2/home.php?action=switch_db). \
1. <geneID>.txt contains the extracted table data containing variants specific to the specified geneID in the Leiden \
Database. Each variant is on its own line and columns are separated by commas. Header labels are included as the first \
line of the file.\
2. <geneID>_MutalizerInput.txt contains only the DNA Change column of <geneID>.txt (one variant per line). This file \
can be directly input to the mutalyzer batch position converter tool by LOVD \
(https://mutalyzer.nl/batchPositionConverter)")
group = parser.add_mutually_exclusive_group()
group.add_argument("-d", "--debug", action="store_true",
                   help="When errors are encountered, a full stack traceback is printed.")
group.add_argument("-g", "--availableGenes", action="store_true", help="A list of all available genes is printed.")
group.add_argument("-a", "--all", action="store_true",
                   help="Extract data for all available genes in the Leiden Database.")
parser.add_argument("geneID", help="Gene ID or multiple geneIDs to retrieve from the Leiden Database.", nargs="*")

args = parser.parse_args()


# User has specified the available genes option, print a list of all available genes.
if args.availableGenes:
    print("---> IN PROGRESS...")
    print("\n".join(LeidenDatabase.get_available_genes()))

# User has specified the all option, so extract data from all genes available on the Leiden Database
elif args.all:
    for gene in LeidenDatabase.get_available_genes():
        try:
            print("---> " + gene + ": IN PROGRESS...")
            LeidenDatabase.save_gene_data(gene)
            print("---> " + gene + ": COMPLETE")
        except:
            LeidenDatabase.print_errors(args, gene)

# The user has not specified all, process their arguments
else:
    # No arguments passed
    if len(args.geneID) == 0:
        print("---> NO GENES PROCESSED: Must pass at least one file or use the --all option")

    # Process each gene ID the user has passed
    else:
        for gene in args.geneID:
            try:
                print("---> " + gene + ": IN PROGRESS...")
                LeidenDatabase.save_gene_data(gene)
                print("---> " + gene + ": COMPLETE")
            except:
                LeidenDatabase.print_errors(args, gene)
