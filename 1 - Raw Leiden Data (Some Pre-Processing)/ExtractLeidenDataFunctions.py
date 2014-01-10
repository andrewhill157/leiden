from LeidenDatabase import *
import traceback
import time
from collections import namedtuple

# TODO update documentation
def get_remapping_results(gene_ids, id_numbers):
    """

    @param gene_ids:
    @type gene_ids:
    @param id_numbers:
    @type id_numbers:
    @return:
    @rtype:
    """

    remapper = VariantRemapper()
    results = []

    for i in range(0, len(gene_ids)):
        print('    ---> ' + gene_ids[i])

        if id_numbers[i] > 0:
            while remapper.entries_remaining_in_batch(id_numbers[i]) > 0:
                time.sleep(0.5)

            results.append(remapper.get_batch_results(id_numbers[i]))
            print('        ---> Complete')

        else:
            results.append([])
            print('        ---> ERROR')

    return results

# TODO update documentation
def format_output_text(header, table_data, remapping):
    """

    @param header:
    @type header:
    @param table_data:
    @type table_data:
    @param remapping:
    @type remapping:
    @return:
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


# TODO update documentation
def write_output_file(gene_id, output_data):
    """

    @param gene_id:
    @type gene_id:
    @param output_data:
    @type output_data:
    @return:
    @rtype:
    """

    # Constants for file delimiters
    file_extension = '.txt'
    row_delimiter = '\n'
    column_delimiter = '\t'

    filename = "".join([gene_id, file_extension])

    # write table data to file in Unicode encoding (some characters are not ASCII encodable)
    with open(filename, 'w', encoding='utf-8') as f:
        file_lines = []

        for row in output_data:
            file_lines.append(column_delimiter.join(row))

        f.write(row_delimiter.join(file_lines))


# TODO update documentation
def process_gene(leiden_database, gene_id):
    """

    @param leiden_database:
    @type leiden_database:
    @param gene_id:
    @type gene_id:
    @return:
    @rtype:
    """
    try:
        # Get header data
        print('    ---> Downloading Variant Data...')
        headers = leiden_database.get_table_headers(gene_id)
        hgvs_mutation_column = Utilities.find_string_index(headers, u'DNA\xa0change')

        print('    ---> Processing Variant Data...')
        # Get table data for variants on given gene
        entries = leiden_database.get_table_data(gene_id)
        variants = [x[hgvs_mutation_column] for x in entries]

        print('    ---> Submitting Remapping...')
        # Remap all variants to genomic coordinates
        remapper = VariantRemapper()
        id_number = remapper.submit_variant_batch(variants)

        if id_number > 0:
            print('    ---> Remapping Submitted.')
        else:
            print('    ---> ERROR: Batch Remapping Failed.')
    except:
        id_number = -1
        entries = []
        headers = []

    process_gene_results = namedtuple('process_gene_results', 'id_number table_data, header_data')
    return process_gene_results(id_number, entries, headers)

# TODO add unit test
def match_column_order(header_template, data_header, table_data):
    """
    Given data and column labels, matches the order of columns in the data to a provided template. header_template
    represents the desired order of the columns and data_header represents the order of columns in provided table_data.
    Columns in table_data are reordered to match the order indicated by header_template.

    @param header_template: desired order of columns in table_data. Must be equal in length to data_header and contain \
    the same items.
    @type header_template: list
    @param header: order of columns in table_data. Must be equal in length to header_template and contain the same items.
    @type header: list
    @param table_data: list of lists. Each inner list represents a row of data and must match the length of \
    data_header and header_template.
    @type: list of lists
    @return: table_data with columns reordered to match the order indicated by the header_template
    @raise: ValueError if contents of lists do not match when sorted
    """

    if(sorted(data_header) == sorted(header_template)):
        for i in range(0, len(header_template)):
            # Find
            header_index = Utilities.find_string_index(data_header, header_template[i])

            # Swap respective rows in table_data if positioning is incorrect
            if i is not header_index:
                table_data = [Utilities.swap(row, i, header_index) for row in table_data]
                data_header = Utilities.swap(data_header, i, header_index)

    else:
        raise ValueError('data_header and header_template do not contain the same elements')

    return table_data


# TODO update documentation
def get_errors(gene_id):
    """

    @param gene_id:
    @return:
    @rtype:
    """

    return '    ---> ERROR: No Entries Found. Please Verify.'
