from LeidenDatabase import *
import time
import copy
from collections import namedtuple

def get_remapping_results(batch_id_numbers):
    """
    Returns batch remapping results for all batch_id_numbers in list. Assumes that all id_numbers in batch_id_numbers
    are valid and are ids for jobs that have already been submitted.

    @param batch_id_numbers: batch remapping id numbers
    @type batch_id_numbers: list
    @return: remapped variants
    @rtype: list
    """

    remapper = VariantRemapper()
    results = []

    for id_number in batch_id_numbers:
        if id_number > 0:
            while remapper.entries_remaining_in_batch(id_number) > 0:
                time.sleep(0.5)

            results.append(remapper.get_batch_results(id_number))

        else:
            results.append([])

    return results


def format_output_text(header, table_data, remapping):
    """
    Combine data from header, table_data, and remapping to form single list of lists.

    @param header: column labels (must match number of columns in table_data)
    @type header: list
    @param table_data: list of lists where inner lists represent rows of data in a table
    @type table_data: list of lists
    @param remapping: genomic mappings of variants (must match number of entries in table_data)
    @type remapping: named tuples of lists. Entries are chromosome_number, coordinate, ref, and alt.
    @return: combined data from headers, table_data, and remapped_variants
    @rtype: list of lists
    """
    # Copy to prevent input variables from being modified
    result = Utilities.deep_copy(table_data)
    result_header = Utilities.deep_copy(header)

    if len(remapping) > 0:
        # Insert new column headers for remapping results
        remapping_start_index = Utilities.find_string_index(header, u'DNA\xa0change') + 1
        result_header[remapping_start_index:remapping_start_index] = ['Chromosome Number', 'Coordinate', 'Ref', 'Alt']

        # Insert remapping result columns
        for i in range(0, len(remapping.chromosome_number)):
            result[i][remapping_start_index:remapping_start_index] = \
                [remapping.chromosome_number[i], remapping.coordinate[i], remapping.ref[i], remapping.alt[i]]

    # Combine all data into one list of lists
    result.insert(0, result_header)
    return result

# TODO clean up extraneous code
def format_vcf_text(header, table_data, remapping):
    """
    Create formatted VCF file data from remapped variants.

    @param header: column labels (must match number of columns in table_data)
    @type header: list
    @param table_data: list of lists where inner lists represent rows of data in a table
    @type table_data: list of lists
    @param remapping: genomic mappings of variants (must match number of entries in table_data)
    @type remapping: named tuples of lists. Entries are chromosome_number, coordinate, ref, and alt.
    @return: combined data from headers, table_data, and remapped_variants
    @rtype: list of lists
    """

    # Initialize with required header information for the VCF input to Variant Effect Predictor
    vcf_text = [
        ['##fileformat=VCFv4.0'],
        ['##INFO=<ID=LAA_CHANGE,Number=1,Type=String,Description="LOVD amino acid change">'],
        ['##INFO=<ID=HGVS_NOTATION,Number=1,Type=String,Description="LOVD HGVS notation describing DNA change">'],
        ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    ]

    # Build remaining rows using remapped variants
    for i in range(0, len(remapping.chromosome_number)):
        # Extract information for tags in VCF file
        amino_acid_change_column = Utilities.find_string_index(header, 'Protein')
        laa_change = table_data[i][amino_acid_change_column]

        hgvs_notation_column = Utilities.find_string_index(header, 'DNA')
        hgvs_notation = table_data[i][hgvs_notation_column]
        try:
            hgvs_notation = hgvs_notation.split(':')[1]
        except:
            print(hgvs_notation)

        vcf_file_row = [remapping.chromosome_number[i], remapping.coordinate[i], '.', remapping.ref[i], remapping.alt[i], '.',
               '.', 'LAA_CHANGE='+laa_change + ';HGVS=' + hgvs_notation]
        vcf_text.append(vcf_file_row)

    return vcf_text


def write_output_file(file_name, output_data):
    """
    Writes output_data to tab-delimited file with one row of data per line.

    @param file_name: name of output file with extension (can include path)
    @type file_name: string
    @param output_data: table data to output to file
    @type output_data: list of lists
    """

    # Constants for file delimiters
    row_delimiter = '\n'
    column_delimiter = '\t'

    # write table data to file in Unicode encoding (some characters are not ASCII encodable)
    with open(file_name, 'w', encoding='utf-8') as f:
        file_lines = []

        for row in output_data:
            file_lines.append(column_delimiter.join(row))

        f.write(row_delimiter.join(file_lines))


def extract_data_and_submit_remap(leiden_database, gene_id):
    """
    Extracts variant table data for given gene in leiden_database and submits variants for remapping to genomic coordinates.

    @param leiden_database: database containing tables of variant data for specified gene_id
    @type leiden_database: LeidenDatabase
    @param gene_id: a string with the Gene ID of the gene to be extracted.
    @type gene_id: string
    @return: namedtuple containing remapping_batch_id_number, table_entries, and column_label entries.
    @rtype: namedtuple
    """
    try:
        print('    ---> Processing Variant Data...')

        # Get header data
        column_labels = leiden_database.get_table_headers()
        hgvs_mutation_column = Utilities.find_string_index(column_labels, u'DNA\xa0change')

        # Get table data for variants on given gene
        table_entries = leiden_database.get_table_data()
        variants = [x[hgvs_mutation_column] for x in table_entries]

        print('    ---> Submitting Remapping...')
        remapper = VariantRemapper()
        remapping_batch_id_number = remapper.submit_variant_batch(variants)

        if remapping_batch_id_number > 0:
            print('    ---> Remapping Submitted.')
        else:
            print('    ---> ERROR: Batch Remapping Failed.')
    except:
        remapping_batch_id_number = -1
        table_entries = []
        column_labels = []

    results = namedtuple('results', 'remapping_batch_id_number table_entries, column_labels')
    return results(remapping_batch_id_number, table_entries, column_labels)