import re
import glob
import traceback


def map_aa_codes(code):
    """
    TODO document
    @param code:
    @rtype: string
    @return:
    """

    # Mapping from three letter amino acid codes to one letter amino acid codes
    one_letter_aa_codes = dict(VAL='V', ILE='I', LEU='L', GLU='E', GLN='Q', ASP='D', ASN='N', HIS='H', TRP='W', PHE='F',
                               TYR='Y', ARG='R', LYS='K', SER='S', THR='T', MET='M', ALA='A', GLY='G', PRO='P', CYS='C',
                               X='*', XAA='*', SCY='*')  # Different notations for stop codons

    # Mapping from one letter amino acid codes to three letter amino acid codes
    three_letter_aa_codes = dict([[v,k] for k, v in one_letter_aa_codes.items() if v != '*'])

    if "*" in code:
        return code

    code_length = len(code)
    if code_length == 3 or code == 'X':
        return one_letter_aa_codes[code]
    elif code_length == 1:
        return three_letter_aa_codes[code]
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
        raw_annotation = remove_parentheses(raw_annotation)
        pattern = re.compile('([A-Za-z\*]{1,3})[\d+]+([A-Za-z\*]{1,3})[A-Za-z\*]*\d*')
        match = re.search(pattern, raw_annotation)

        if match is not None:
            return [match.group(1), match.group(2)]
        else:
            return ['-', '-']

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

def get_tagged_files(tag):
    """
    Returns a list of all .vcf files in the directory containing the specified substring, tag in the filename.
    @param tag: a string which is present in the files you want to include in the returned list
    @rtype: list of strings
    @return: List of all files with a .vcf extension that tag as a substring in the filename.
    """
    return [x for x in glob.glob("*.vcf") if tag in x]

def flag_annotation_file(annotation_file):
    """
    TODO document
    @param annotation_file:
    """
    flag_index = annotation_file.find("_")

    # <original file up to the start of _FINAL.vcf>_FLAGGED.vcf
    results_file = annotation_file[0:flag_index] + "_FLAGGED.vcf"

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
                        print('LOVD: ' + str(laa_change) +'; VEP: ' + str(aa_change))
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

def print_errors(commandline_args, file_name):
    """
    Given the list of command line arguments passed to the script and a file_name, print error messages to the console.

    @param commandline_args: a parser.parse_args() object from the argparse library. Assumed to contain an argument \
    called debug to indicate verbosity of error messages. args.debug is true, full stack traces are printed for \
    errors. A simple error message is printed otherwise.
    @param file_name: a string with the file_name of the file that generated the error during processing.
    """
    if commandline_args.debug:
        print("---> " + file_name + ": ERROR - NOT PROCESSED. STACK TRACE: ")
        tb = traceback.format_exc()
        print(tb)
    else:
        print("---> " + file_name + ": ERROR - NOT PROCESSED. Use --debug option for more details.")

