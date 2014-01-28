from LeidenDatabase import *
import traceback
import time
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
    @param remapped_variants: genomic mappings of variants (must match number of entries in table_data)
    @type remapped_variants: list
    @return: combined data from headers, table_data, and remapped_variants
    @rtype: list of lists
    """

    if len(remapping) > 0:
        # Insert new column headers for remapping results
        remapping_start_index = Utilities.find_string_index(header, u'DNA\xa0change') + 1
        header[remapping_start_index:remapping_start_index] = ['Chromosome Number', 'Coordinate', 'Ref', 'Alt']

        # Insert remapping result columns
        for i in range(0, len(remapping.chromosome_number)):
            table_data[i][remapping_start_index:remapping_start_index] = \
                [remapping.chromosome_number[i], remapping.coordinate[i], remapping.ref[i], remapping.alt[i]]

    # Combine all data into one list of lists
    table_data.insert(0, header)
    return table_data


def write_output_file(file_name, output_data):
    """
    Writes output_data to tab-delimited file with one row of data per line.

    @param file_name: name of output file without extension (can include path)
    @type file_name: string
    @param output_data: table data to output to file
    @type output_data: list of lists
    """

    # Constants for file delimiters
    file_extension = '.txt'
    row_delimiter = '\n'
    column_delimiter = '\t'

    filename = "".join([file_name, file_extension])

    # write table data to file in Unicode encoding (some characters are not ASCII encodable)
    with open(filename, 'w', encoding='utf-8') as f:
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


def match_column_order(header_template, data_header, table_data):
    """
    Given data and column labels, matches the order of columns in the data to a provided template. header_template
    represents the desired order of the columns and data_header represents the order of columns in provided table_data.
    Columns in table_data are reordered to match the order indicated by header_template.

    @param header_template: desired order of columns in table_data. Must be equal in length to data_header and contain
    the same items.
    @type header_template: list
    @param data_header: order of columns in table_data.
    @param table_data: list of lists. Each inner list represents a row of data and must match the length of
    data_header and header_template.
    @type: list of lists
    @return: table_data with columns reordered to match the order indicated by the header_template
    @raise: ValueError if contents of lists do not match when sorted
    """

    if(sorted(data_header) == sorted(header_template)):
        for i in range(0, len(header_template)):
            # Find index of current string from template in data_header
            header_index = Utilities.find_string_index(data_header, header_template[i])

            # Swap respective rows in table_data if positioning is incorrect
            if i is not header_index:
                table_data = [Utilities.swap(row, i, header_index) for row in table_data]
                data_header = Utilities.swap(data_header, i, header_index)

    else:
        raise ValueError('data_header and header_template do not contain the same elements')

    return table_data