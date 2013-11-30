import argparse
from MutalyzerProcessing import *


parser = argparse.ArgumentParser(description="Combines output from the Mutalyzer Batch Position Converter Tool and the \
original extracted data from the Leiden Database entry for a given gene. This produces a file in the called \
<geneID>_MAPPED, which contains the combined information of <geneID>.txt (original extracted Leiden Database data) and \
<geneID>_MutalyzerOutput (contains the mutalyzer output from all variants from the original Leiden Database. These two \
files are required and must be present in the same directory to run the script and must be named in this manner. \
It is assumed that there there are no other .txt files in the directory. All non .txt files are ignored and .txt files\
with an '_' in the name (such as the file output by this script) are ignored. \
Note that the output file is also produced in this same directory. It is assumed that all variants from the original \
data are present in the both files.")
group = parser.add_mutually_exclusive_group()
group.add_argument("-d", "--debug", action="store_true", help="When errors are encountered, a full stack traceback is \
printed.")

group2 = parser.add_mutually_exclusive_group()
group2.add_argument("-a", "--all", action="store_true", help="All files with a .txt extension that do not contain \
MutalyzerOutput or _MAPPED in the filename are processed.")
parser.add_argument("fileNames", help="File name or multiple file names to compile. This should be a text file with \
the original extracted Leiden Database data (the file ACTA1.txt when calling python extractLeidenData.py ACTA1, \
for example)", nargs="*")

args = parser.parse_args()

# User has specified the all option, process all files in the directory
if args.all:
    for files in list_file_in_directory():
        try:
            annotation_input = get_annotation_input(files)
            write_table_to_file(annotation_input, remove_file_extension(files) + "_MAPPED.txt")
            print("---> " + files + ": COMPLETE")
        except:
            print_errors(args, files)

# The user has not specified all, process their arguments
else:
    # No arguments passed
    if len(args.fileNames) == 0:
        print("---> NO FILES PROCESSED: Must pass at least one file or use the --all option")
    # Process each file the user passes as an argument
    else:
        for files in args.fileNames:
            if ".txt" in files and not "_" in files:
                try:
                    annotation_input = get_annotation_input(files, remove_file_extension(files) + "_MutalyzerOutput.txt")
                    write_table_to_file(annotation_input, remove_file_extension(files) + "_MAPPED.txt")
                    print("---> " + files + ": COMPLETE")
                except:
                    print_errors(args, files)
            else:
                print("---> " + files + ": ERROR - NOT PROCESSED. .txt extension must be included with file arguments.")