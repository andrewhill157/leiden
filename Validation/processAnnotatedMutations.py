import argparse
from AnnotationProcessing import *


parser = argparse.ArgumentParser(description="Given the VCF-like file output by annotation of variants, add a column \
    of flags as the first column in the file to summarize information contained in the INFO column. \
    The script will output files in the same directory as the script called <original_file_name>_FLAGGED.vcf \
    For example, the command python processAnnotatedMutations ACTA1_FINAL.vcf would produce a file called \
    ACTA1_FLAGGED.vcf, etc. \
    \
    The output file will contain all the same information as the original file, except the column will now contain a \
    semi-colon delimited list of tags and corresponding values for each mutation: TAG1=VALUE;TAG2=VALUE;TAG3=VALUE \
    Note that the original columns in the file will all be shifted right one column in the output.\
    Tag and value definitions in resulting output file: \
    CONCORDANCE= \
        - CONCORDANT: The amino acid changes predicted by LAA_CHANGE and AA_CHANGE in the input file match \
        - NOT_CONCORDANT: The amino acid changes predicted by LAA_CHANGE and AA_CHANGE in the input file do not match \
        - ERROR: There was an error processing the LAA_CHANGE and/or AA_CHANGE entry \
        - NO_AA_CHANGE: No amino acid change is predicted by the annotation \
    \
    SEVERE_IMPACT= \
        - The value of SEVERE_IMPACT in the original file is simply copied as the value of this tag for convenience \
    \
    HGMD= \
        - OVERLAPS - An HGMD_MUT entry is present in the original file \
        - NOT_IN_HGMD - No HGMD_MUT entry is present in the original file \
    \
    ALLELE_FREQUENCY=  \
        - LOW_FREQUENCY: Overall allele frequency,(AC_MAC26K/AN_MAC26K)*100, is less than or equal to 0.5% \
        - HIGH_FREQUENCY: Overall allele frequency,(AC_MAC26K/AN_MAC26K)*100, is greater than 0.5% \
        - NOT_IN_26K_DATABASE: There were no AC_MAC26K or AN_MAC26K entries for this mutation")

group = parser.add_mutually_exclusive_group()
group.add_argument("-d", "--debug", action="store_true",
                   help="When errors are encountered, a full stack traceback is printed.")
group.add_argument("-a", "--all", action="store_true",
                   help="Process all vcf files in the given directory. Only files containing FINAL in the file name \
                   will be processed.")
parser.add_argument("file_name", help="File or files (including extension) to process.", nargs="*")

args = parser.parse_args()

if args.all:
    for files in get_tagged_files("FINAL"):
        try:
            flag_annotation_file(files)
            print("---> " + files + ": COMPLETE")
        except:
            print_errors(args, files)
else:
    for files in args.file_name:
        try:
            flag_annotation_file(files)
            print("---> " + files + ": COMPLETE")
        except: 
            print_errors(args, files)
    


            