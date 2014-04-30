"""
Andrew Hill
MacArthur Lab - 2014

Functions for handling file data input and output.
"""


def format_vcf_text(vcf_format_variants, info_column_entries):
    """
    Create formatted VCF file data from remapped variants.

    @param vcf_format_variants: list of tuples that contain chromosome_number, coordinate, ref, alt in that order
    @type: list of tuples of strings
    @param info_column_entries: dictionary where keys are tags to place in info column. The values are tuples that
    contain data_type, description, and a list of values (one per entry in vcf_format_variant)
    @type info_column_entries: dict
    @return: list of lists where each inner list is a row of the VCF file, and each element in the inner list represents
    the value of the respective column in a given row.
    @rtype: list of lists
    """
    # Initialize with required header information for the VCF input to Variant Effect Predictor
    vcf_text = [
        ['##fileformat=VCFv4.0'],
        ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    ]

    # Insert INFO tag descriptions into VCF header
    tag_keys = info_column_entries.keys()

    for tag in tag_keys:
        tag_header_text = ''.join(['##INFO=<ID=', tag, ',Number=1,Type=', info_column_entries[tag][0], ',Description="', info_column_entries[tag][1], '">'])
        vcf_text.insert(1, [tag_header_text])

    # Format VCF body text
    for i in range(0, len(vcf_format_variants)):

        chromosome_number, coordinate, ref, alt = vcf_format_variants[i]

        tag_list = []
        for tag in tag_keys:
            tag_value = info_column_entries[tag][2][i].replace(' ', '')
            tag_list.append(tag + '=' + tag_value)

        vcf_file_row = [chromosome_number,
                        coordinate,
                        '.',
                        ref,
                        alt, '.',
                        '.',
                        ';'.join(tag_list)]

        vcf_text.append(vcf_file_row)

    return vcf_text


def read_table_from_file(file_name, column_delimiter='\t'):
    """
    Returns table of data from tab-delimited file. Alternate delimiter can be specified.

    @param file_name: path to tab-delimited file wit
    @type file_name: string
    @param column_delimiter: column delimiter (tab by default)
    @type column_delimiter: string
    @return:table of data from tab-delimited file
    @rtype: 2D list of strings (1st dimension is rows, 2nd is columns)
    """

    with open(file_name, 'r') as f:
        file_text = f.read().splitlines()

    return [line.split('\t') for line in file_text]


def write_table_to_file(file_name, table, column_delimiter='\t'):
    """
    Writes table to tab-delimited file. Alternate delimiter can be specified.

    @param file_name: name of output file with extension (can include path)
    @type file_name: string
    @param column_delimiter: column delimiter (tab by default)
    @type column_delimiter: string
    @param table: table data to output to file
    @type table: 2D list of strings (1st dimension is rows, 2nd is columns)
    """

    row_delimiter = '\n'

    file_lines = []

    for row in table:
        file_lines.append(column_delimiter.join(row))

    output_text = row_delimiter.join(file_lines)

    # Ensure unicode strings are encoded before writing
    if isinstance(output_text, unicode):
        output_text = output_text.encode('utf-8')

    with open(file_name, 'w') as f:
        f.write(output_text)