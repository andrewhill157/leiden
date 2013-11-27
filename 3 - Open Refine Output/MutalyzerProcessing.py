import re
import os
import traceback


def get_combined_data(raw_leiden_file, mutalyzer_output_file):
    """
    Returns the combined data from the raw leiden data output file and the remapped mutalyzer output file for a given
    gene.
    @param raw_leiden_file: a text file containing the raw extracted Leiden Database data (as output by the \
    extractLeidenData.py script. For example, ACTA1.txt as output by python extractLeidenData.py ACTA1.
    @param mutalyzer_output_file: Mutalyzer batch position converter output text file from variant remapping \
    (as obtained when submitting variants to Mutalyzer). For example, the resulting output file when submitting \
    ACTA1_MutalyzerInput.txt produced by python extractLeidenData.py ACTA1.
    @rtype: list of lists of strings
    @return: list of lists where each inner list is a row of the data combined from the raw_leiden_file and \
    mutalyzer_output_file. The order of the columns matches the order of columns in raw_leiden_file from left to right.\
    The two columns of data from mutalyzer_output file (Errors and Chromosomal Variant) are inserted at the as the \
    second and third columns of the combined data respectively. The first inner list contains column labels for all \
    columns in the combined data. The column labels from mutalyzer output are inserted in the same positions as their \
    respective columns of data.
    """
    with open(raw_leiden_file) as leidenData, open(mutalyzer_output_file) as mutalyzerOutput:
        # Files are tab delimited
        column_delimiter = "\t"

        # Identify the Errors and Chromosomal Variant columns in the Mutalyzer Output file
        column_labels = mutalyzerOutput.readline().split(column_delimiter)

        error_column = find_string_index(column_labels, "Errors")
        variant_column = find_string_index(column_labels, "Chromosomal Variant")

        errors = [column_labels[error_column]]
        mappings = [column_labels[variant_column]]

        # Store all entries in from the mutalyzer output
        mutalyzer_line_count = 0
        for line in mutalyzerOutput:
            mutalyzer_line_count += 1  # track number of lines in file

            columns = line.split(column_delimiter)
            errors.append(columns[error_column])
            mappings.append(columns[variant_column])

        # Insert the mutalyzer output into the data
        combined_data = []
        i = 0
        for line in leidenData:
            # Error checking for equal number of variants
            if i > mutalyzer_line_count:
                raise ValueError("Mutalyzer output file and Leiden Data file contain a different number of variants!")
            temp = line.split(column_delimiter)
            temp.insert(2, errors[i])
            temp.insert(3, mappings[i])
            combined_data.append(temp)
            i += 1
            
        # Error checking for equal number of variants
        if i <= mutalyzer_line_count:
            raise ValueError("Mutalyzer output file and Leiden Data file contain a different number of variants!")

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
        column_delimiter = "\t"

        for lists in table:
            f.write(column_delimiter.join(lists))


def remove_file_extension(file_name):
    """
    Removes the file extension from a given file.
    @param file_name: name of the file whose extension is to be removed. Files without extensions returned unmodified.
    @rtype: string
    @return: file_name with no file extension ("." character is also removed). For example an input of myFile.txt \
    would return a new string, myFile
    """
    return os.path.splitext(file_name)[0]


def get_annotation_input(raw_leiden_data_file, mutalyzer_output_file):
    """
    Given a raw leiden database file and a mutalyzer output file, returns the data with remapped variants inserted
    as new columns in the data in a VCF-like format (ERRORS, CHROMOSOME NUMBER, REF, and ALT are inserted columns).
    Assumes that the two files have the same number of variants listed.
    @param raw_leiden_data_file: path to .txt file with the raw leiden data output file (ACTA1.txt as output by python \
    extractLeidenData.py ACTA1, for example).
    @param mutalyzer_output_file: path to .txt file with the output from the mutalyzer remapping tool (output from \
    Mutalyzer when submitting ACTA1_MutalyzerInput.txt as produced by python extractLeidenData.py ACTA1, for example).
    @rtype: list of lists of strings
    @return: combined data from the two files with the remapped variants from Mutalyzer output files inserted under \
    the column label Chromosomal Variant in addition to the 4 VCF format columns for given variants (ERRORS, \
    CHROMOSOME NUMBER, REF, and ALT). All information from the raw_leiden_data_file is preserved. Inner lists of \
    returned value are rows of the data and indices in the inner lists are column entries from left to right (in order \
    of increasing index number).
    """

    # Do not want to process files with an _ in name (avoids trying to process _MAPPED and _MutalyzerOutputFiles
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


