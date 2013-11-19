import re
import os
import argparse
import traceback


def get_combined_data(raw_leiden_file, mutalyzer_output_file):
    """
    Returns the combined data from the raw leiden data output file and the remapped mutalyzer output file for a given
    gene.
    @param raw_leiden_file: a text file containing the raw extracted Leiden Database data (as output by the \
    extractLeidenData.py script. For example, ACTA1.txt as output by python extractLeidenData.py ACTA1.
    @param mutalyzer_output_file: Mutalyzer batch position converter output text file from variant remapping \
    (as obtained when submitting variants to Mutalyzer). For example, the resulting output file when submitting \
    ACTA1_MutalizerInput.txt produced by python extractLeidenData.py ACTA1.
    @rtype: list of lists of strings
    @return: list of lists where each inner list is a row of the data combined from the raw_leiden_file and \
    mutalyzer_output_file. The order of the columns matches the order of columns in raw_leiden_file from left to right.\
    The two columns of data from mutalyzer_output file (Errors and Chromosomal Variant) are inserted at the as the \
    second and third columns of the combined data respectively. The first inner list contains column labels for all \
    columns in the combined data. The column labels from mutalyzer output are inserted in the same positions as their \
    respective columns of data.
    """
    with open(raw_leiden_file) as leidenData, open(mutalyzer_output_file) as mutalyzerOutput:

        # Identify the Errors and Chromosomal Variant columns in the Mutalyzer Output file
        column_labels = mutalyzerOutput.readline().split("\t")
        error_column = find_string_index(column_labels, "Errors")
        variant_column = find_string_index(column_labels, "Chromosomal Variant")

        errors = [column_labels[error_column]]
        mappings = [column_labels[variant_column]]

        # Store all entries in from the mutalyzer output
        for line in mutalyzerOutput:
            columns = line.split("\t")
            errors.append(columns[error_column])
            mappings.append(columns[variant_column])

        # Insert the mutalyzer output into the data
        combined_data = []
        i = 0
        for line in leidenData:
            temp = line.split(",")

            temp.insert(2, errors[i])
            temp.insert(3, mappings[i])
            combined_data.append(temp)
            i += 1

        return combined_data


def find_string_index(string_list, search_string):
    """
    Given a list of strings and a string to search for, returns the first index of the first element in the list that
    contains the search string. Note that the comparison is case sensitive and the element at a given index must only
    contain the search element, not exactly match the search element.
    @param string_list: list of strings
    @param search_string: a string to search for in elements of string_list
    @rtype: number
    @return: index of the first instance of the search string in an element of the list. Returns -1 if the search \
    string is not found.
    """
    i = 0
    for entry in string_list:
        if search_string in entry:
            return i
        else:
            i += 1
    return -1


def get_chromosome_number(mapping):
    """
    Given a mapping of the form NC_######.#:g.<HGVS notation>, where NC_###### is the chromosome reference ID,
    returns tha chromosome number from the mutation mapping. For example, calling on NC_000001.2:g.524264A>C 
    would return 1. 
    The mapping is assumed to be in valid HGVS notation.
    @param mapping: string with the mapping of a mutation in HGVS notation as shown above.
    @rtype: string
    @return: chromosome number of the mutation as a string
    """
    m = re.search('([0]+)([1-9][0-9]?)([.])', mapping)
    return m.group(2)


def get_coordinates(mapping):
    """
    Given a mapping of the form NC_######.#:g.<HGVS notation>, where NC_###### is the chromosome reference ID,
    returns tha coordinate from the mutation mapping. For example, calling on NC_0000001.2:g.524264A>C 
    would return 524264. 
    The mapping is assumed to be in valid HGVS notation.
    @param mapping: string with the mapping of a mutation in HGVS notation as shown above.
    @rtype: string
    @return: coordinates of the mutation as a string
    """
    m = re.search('([g][.])([0-9]+[_]?[0-9]+)', mapping)
    return m.group(2)


def get_ref(mapping):
    """
    Given a mapping of the form NC_######.#:g.<HGVS notation>, where NC_###### is the chromosome reference ID,
    returns tha REF (reference) base from the mutation mapping. For example, calling on NC_0000001:g.524264A>C 
    would return A. 
    The mapping is assumed to be in valid HGVS notation.
    @param mapping: string with the mapping of a mutation in HGVS notation as shown above.
    @rtype: string
    @return: reference base of the mutation (base before the mutation) as a string
    """
    m = re.search('([A-Z])([>])([A-Z])', mapping)
    return m.group(1)


def get_alt(mapping):
    """
    Given a mapping of the form NC_######.#:g.<HGVS notation>, where NC_###### is the chromosome reference ID,
    returns tha ALT (reference) base from the mutation mapping. For example, calling on NC_0000001:g.524264A>C 
    would return C. 
    The mapping is assumed to be in valid HGVS notation.
    @param mapping: string with the mapping of a mutation in HGVS notation as shown above.
    @rtype: string
    @return: alternate base of the mutation (base after the mutation) as a string
    """
    m = re.search('([A-Z])([>])([A-Z])', mapping)
    return m.group(3)


def convert_to_vcf(combined_data):
    """
    Given the combined Leiden Data and Mutalyzer Output (remapped variants), convert the Chromosomal Variant column
    to VCF format (Chromosome Number, Coordinate, Ref, Alt).
    @param combined_data: combined Leiden and Mutalyzer output data as a list of lists, where each inner list contains \
    data for a given row of the combined data table. Columns at the entries ordered left to right. The data must \
    contain the column labels as the first inner list and it must contain a column named Chromosomal Variant.
    @rtype: list of lists strings
    @return: A new list with the Chromosomal Variant column split into four seperate columns. The first of these \
    columns is in the original index of the Chromosomal Variant column in combined_data and the four inserted columns \
    are CHROMOSOME NUMBER, COORDINATE, REF, and ALT from lowest to highest index. Currently, only SNPs are being \
    processed, all other variant types are left unchanged. Variants with no remapped notation due to errors in \
    remapping are returned unchanged. For unchanged entries, the inserted columns will contain an empty string.
    """
    mapping_index = find_string_index(combined_data[0], "Chromosomal Variant")
    vcf_mapping = []
    combined_data[0].insert(mapping_index, "ALT")
    combined_data[0].insert(mapping_index, "REF")
    combined_data[0].insert(mapping_index, "COORDINATE")
    combined_data[0].insert(mapping_index, "CHROMOSOME NUMBER")
    vcf_mapping.append(combined_data[0])
    for row in combined_data[1:]:
        mapping = row[mapping_index]

        if mapping is not "":
            chromosome = get_chromosome_number(mapping)

            if ">" in mapping:
                coordinates = get_coordinates(mapping)
                ref = get_ref(mapping)
                alt = get_alt(mapping)

            else:
                # Not yet processing anything more complicated than SNPs
                coordinates = mapping
                ref = ""
                alt = ""
        else:
            chromosome = ""
            coordinates = ""
            ref = ""
            alt = ""

        row.insert(mapping_index, alt)
        row.insert(mapping_index, ref)
        row.insert(mapping_index, coordinates)
        row.insert(mapping_index, chromosome)
        vcf_mapping.append(row)

    return vcf_mapping


def write_table_to_file(table, file_name):
    """
    Writes a list of lists to a file with a specified file name. 
    Elements of table (lists themselves) are each written to their own line in the file with a tab character separating
    each entry.
    @param table: list whose elements are lists of data elements such as strings, ints, etc.
    @param file_name: File name of the file to write table to. Must be a valid file name with extension.
    """
    with open(file_name, 'w') as f:
        for lists in table:
            f.write("\t".join(lists))


def remove_file_extension(file_name):
    """
    Removes the file extension from a given file. 
    @param file_name: name of the file whose extension is to be removed. Files without extensions returned unmodified.
    @rtype: string
    @return: file_name with no file extension ("." character is also removed). For example an input of myFile.txt \
    would return a new string, myFile
    """
    return os.path.splitext(file_name)[0]


def get_annotation_input(raw_leiden_data_file):
    """
    TODO document
    @param raw_leiden_data_file:
    """

    # Do not want to process files with an _ in name (avoids trying to process _MAPPED and _MutalizerOutputFiles
    if ".txt" in raw_leiden_data_file and "_" not in raw_leiden_data_file:
        mutalyzer_output_file = remove_file_extension(raw_leiden_data_file) + "_MutalizerOutput.txt"
        combined_data = get_combined_data(raw_leiden_data_file, mutalyzer_output_file)

        return convert_to_vcf(combined_data)


def list_file_in_directory():
    """
    Returns a list of all files in the current directory (the directory the script is running in).
    @rtype: list of strings
    @return: list of all files in the current directory (the directory the script is running in).
    """
    return os.listdir(os.path.dirname(os.path.abspath(__file__)))


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


"""
Command line tool for extracting data from the leiden database.
Execute script with -h option for details on arguments and description. 
"""
parser = argparse.ArgumentParser(description="Combines output from the Mutalyzer Batch Position Converter Tool and the \
original extracted data from the Leiden Database entry for a given gene. This produces a file in the called \
<geneID>_MAPPED, which contains the combined information of <geneID>.txt (original extracted Leiden Database data) and \
<geneID>_MutalizerOutput (contains the mutalizer output from all variants from the original Leiden Database. These two \
files are required and must be present in the same directory to run the script and must be named in this manner. \
Note that the output file is also produced in this same directory. It is assumed that all variants from the original \
data are present in the both files.")
group = parser.add_mutually_exclusive_group()
group.add_argument("-d", "--debug", action="store_true", help="When errors are encountered, a full stack traceback is \
printed.")
group.add_argument("-a", "--all", action="store_true", help="All files with a .txt extension that do not contain \
MutalizerOutput or _MAPPED in the filename are processed.")
parser.add_argument("fileNames", help="File name or multiple file names to compile. This should be a text file with \
the original extracted Leiden Database data (the file ACTA1.txt when calling python extractLeidenData.py ACTA1, \
for example)", nargs="*")

args = parser.parse_args()

# User has specified the aal option, process all files in the directory
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
            if ".txt" in files:
                try:
                    annotation_input = get_annotation_input(files)
                    write_table_to_file(annotation_input, remove_file_extension(files) + "_MAPPED.txt")
                    print("---> " + files + ": COMPLETE")
                except:
                    print_errors(args, files)
            else:
                print("---> " + files + ": ERROR - NOT PROCESSED. .txt extension must be included with file arguments.")