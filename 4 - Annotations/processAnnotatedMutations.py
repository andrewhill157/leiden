import re
import sys
import glob
import argparse

def map_aa_codes(code):
    """
    TODO document
    @param code:
    @rtype: string
    @return:
    """

    # Mapping from three letter amino acid codes to one letter amino acid codes 
    one_letter_aa_codes = dict(VAL='V', ILE='I', LEU='L', GLU='E', GLN='Q', ASP='D', ASN='N', HIS='H', TRP='W', PHE='F',
                               TYR='Y', ARG='R', LYS='K', SER='S', THR='T', MET='M', ALA='A', GLY='G', PRO='P', CYS='C')
    
    # Mapping from one letter amino acid codes to three letter amino acid codes 
    three_letter_aa_codes = dict([[v,k] for k, v in one_letter_aa_codes.items()])

    if "*" in code:
        return code
    code_length = len(code)
    if code_length is 1:
        return three_letter_aa_codes[code]
    elif code_length is 3:
        return one_letter_aa_codes[code]
    else:
        return "ERROR"

def get_annotation_info(line):
    """
    TODO document
    @param line:
    @rtype: list of strings
    @return:
    """

    # Lines in original file are tab separated and info is in 8th column
    info_column = line.split("\t")[7]
    return info_column.split(";")


def get_tagged_entry(annotation_info, tag):
    """
    TODO document
    @param annotation_info:
    @param tag:
    @rtype: string
    @return:
    """

    for entry in annotation_info:
        if entry.startswith(tag):
            return remove_tag(entry, tag)
    return ""


def remove_tag(info_text, tag):
    """
    TODO document
    @param: info_text
    @param: tag
    @rtype: string
    @return:
    """

    tag_length = len(tag)
    return info_text[tag_length:].strip()


def remove_parentheses(annotation_text):
    """
    TODO document
    @param annotation_text:
    @rtype: string
    @return:
    """

    opening_parenthesis_index = annotation_text.find("(")
    pdot_notation_index = annotation_text.find(".")
    
    if opening_parenthesis_index > 0:
        closing_parenthesis_index = annotation_text.find(")")
        return annotation_text[opening_parenthesis_index+1:closing_parenthesis_index]
    elif pdot_notation_index > 0:
        return annotation_text[pdot_notation_index+1:]
    else: return annotation_text

        
def get_laa_change(annotation_info):
    """
    TODO document
    @param annotation_info:
    @rtype: list of strings
    @return:
    """

    laa_change_tag = "LAA_CHANGE="
    raw_annotation = get_tagged_entry(annotation_info, laa_change_tag)
    
    if "-" in raw_annotation or "?" in raw_annotation or "=" in raw_annotation:
        return ""  # no change or unknown change in LAA
    else:
        raw_annotation = remove_parentheses(raw_annotation);
        return re.split("[0-9]+", raw_annotation)
                

def get_aa_change(annotation_info):
    """
    TODO document
    @param annotation_info:
    @rtype: list of strings
    @return:
    """

    aa_change_tag = "AA_CHANGE="
    raw_annotation = get_tagged_entry(annotation_info, aa_change_tag)
    raw_annotation = raw_annotation.split(",")

    for change in raw_annotation:
        if "/" in change:
            return change.split("/")
    return ""


def is_concordant(laa_change, aa_change):
    """
    TODO document
    @param laa_change:
    @param aa_change:
    @rtype: string
    @return:
    """
    try:
        if len(laa_change) is not 2 or len(aa_change) is not 2:
            return "ERROR"
        laa_change = [x.upper() for x in laa_change]
        aa_change = [x.upper() for x in aa_change]
        laa_change = [map_aa_codes(x) for x in laa_change]
        if laa_change == aa_change:
            return "OK"
        else: 
            return "FAIL"
    except: 
        return "ERROR"


def get_severe_impact(annotation_info):
    """
    TODO document
    @param annotation_info:
    @rtype: string
    @return
    """
    severe_impact_flag = "SEVERE_IMPACT="
    return get_tagged_entry(annotation_info, severe_impact_flag)


def get_allele_frequency(annotation_info):
    """
    TODO document
    @param annotation_info:
    @rtype:
    @return:
    """

    ac_flag = "AC_MAC26K="
    an_flag = "AN_MAC26K="
    
    ac = get_tagged_entry(annotation_info, ac_flag)
    an = get_tagged_entry(annotation_info, an_flag)
    
    if ac == "" or an == "":
        return ""
    else:
        ac = float(ac)
        an = float(an)
        return (ac/an)*100


def flag_annotation_file(annotation_file):
    """
    TODO document
    @param annotation_file:
    """

    # Do not process any flagged files in the directory
    if "FLAGGED_" not in annotation_file:
        results_file = "".join(["FLAGGED_", annotation_file])

        with open(annotation_file) as entries, open(results_file, 'w') as results:
            for line in entries:
                # Only process flags for lines that are not header information
                if "#" not in line:
                    flags = []
                    # STEP 1 - Check that LAA_CHANGE= and AA_CHANGE= entries are the same (that entries are concordant)
                    annotation_info = get_annotation_info(line)
                    laa_change = get_laa_change(annotation_info)
                    aa_change = get_aa_change(annotation_info)
                    
                    concordance_flag = "CONCORDANCE="
                    if laa_change is "" or aa_change is "":
                        flags.append("".join([concordance_flag, "NO_AA_CHANGE"]))
                        
                    else:
                        concordance = is_concordant(laa_change, aa_change)
                        if concordance == "OK":
                            flags.append("".join([concordance_flag, "CONCORDANT"]))
                        elif concordance == "FAIL":
                            flags.append("".join([concordance_flag, "NOT_CONCORDANT"]))
                        elif concordance == "ERROR":
                            flags.append("".join([concordance_flag, "ERROR"]))
                    
                    # STEP 2 - Extract the SEVERE_IMPACT entry
                    severe_impact_flag = "SEVERE_IMPACT="
                    severe_impact = get_severe_impact(annotation_info)
                    flags.append("".join([severe_impact_flag, severe_impact]))
                        
                    # STEP 3 - Check for overlap with HGMD
                    hgmd_flag = "HGMD="
                    if "HGMD_MUT" not in line:
                        flags.append("".join([hgmd_flag, "NOT_IN_HGMD"]))
                    else: 
                        flags.append("".join([hgmd_flag, "OVERLAPS"]))
                            
                    # Step 4 - Check for allele frequency information (above or below overall 0.5% frequency
                    allele_frequency_flag = "ALLELE_FREQUENCY="
                    allele_frequency = get_allele_frequency(annotation_info)
                    if allele_frequency == "":
                        flags.append("".join([allele_frequency_flag, "NOT_IN_26K_DATABASE"]))
                    elif allele_frequency > 0.5:
                        flags.append("".join([allele_frequency_flag, "HIGH_FREQUENCY"]))
                    else: 
                        flags.append("".join([allele_frequency_flag, "LOW_FREQUENCY"]))
                        
                    # Write out lines to new file with [LIST_OF_FLAGS] inserted as the new first column
                    results.write(";".join(flags))
                    results.write("\t")
                    results.write(line)
                
                # Copy header data unchanged from original file
                elif line.startswith("##"):
                    results.write(line)
                    
                # Copy column labels with new column, FLAGS, inserted as label for first column
                elif line.startswith("#"):
                    results.write("".join([line[0], "FLAGS", "\t", line[1:]]))
        
"""
BEGIN SCRIPT
"""

parser = argparse.ArgumentParser(description="Given the VCF-like file output by annotation of variants, add a column \
    of flags as the first column in the file to summarize information contained in the INFO column. \
    The script will output files in the same directory as the script called <original_file_name>_FLAGGED.vcf \n\
For example, the command shown in number one above would produce files called FLAGGED_myAnnotatedMutations.vcf, etc. \
\n\
The output file will contain all the same information as the original file, except the column will now contain a \
semi-colon delimited list of tags and corresponding values for each mutation: TAG1=VALUE;TAG2=VALUE;TAG3=VALUE \
Note that the original columns in the file will all be shifted right one column in the output.\
Tag and value definitions in resulting output file: \n\
CONCORDANCE= \n\
    - CONCORDANT: The amino acid changes predicted by LAA_CHANGE and AA_CHANGE in the input file match \n\
    - NOT_CONCORDANT: The amino acid changes predicted by LAA_CHANGE and AA_CHANGE in the input file do not match \n\
    - ERROR: There was an error processing the LAA_CHANGE and/or AA_CHANGE entry \n\
    - NO_AA_CHANGE: No amino acid change is predicted by the annotation \n\
\n\
SEVERE_IMPACT= \n\
    - The value of SEVERE_IMPACT in the original file is simply copied as the value of this tag for convenience \n\
\n\
HGMD= \n\
    - OVERLAPS - An HGMD_MUT entry is present in the original file \n\
    - NOT_IN_HGMD - No HGMD_MUT entry is present in the original file \n\
\n\
ALLELE_FREQUENCY=  \n\
    - LOW_FREQUENCY: Overall allele frequency,(AC_MAC26K/AN_MAC26K)*100, is less than or equal to 0.5% \n\
    - HIGH_FREQUENCY: Overall allele frequency,(AC_MAC26K/AN_MAC26K)*100, is greater than 0.5% \n\
    - NOT_IN_26K_DATABASE: There were no AC_MAC26K or AN_MAC26K entries for this mutation")

group = parser.add_mutually_exclusive_group()
group.add_argument("-d", "--debug", action="store_true",
                   help="When errors are encountered, a full stack traceback is printed.")
group.add_argument("-a", "--all", action="store_true",
                   help="Process all vcf files in the given directory. There can be no other types of VCF files in the \
                   directory except those containing FLAGGED in the name, which are ignored.")
parser.add_argument("file_name", help="File or files (including extension) to process.", nargs="*")

args = parser.parse_args()

arguments = len(sys.argv)
if arguments > 1:
    for files in sys.argv[1:]:
        try:
            flag_annotation_file(files)
        except: 
            print("Error Processing File: " + files)
else:
    for files in glob.glob("*.vcf"):
        try:
            flag_annotation_file(files)
        except: 
            print("Error Processing File: " + files)
    


            