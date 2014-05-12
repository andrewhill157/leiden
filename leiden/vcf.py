import re
import pandas as pd

VCF_HEADER_PREFIX = '#'
VCF_DELIMTER = '\t'
FORMAT_DELIMITER = '|'


def get_vcf_dict_from_file(vcf_file):
    """
    Returns a list, where each each element is a nested dictionary providing access to progressively lower-level
    information from a VCF file.

    The first level is the VCF column name (all caps). The info column is a nested dictionary where tag names can be
    accessed. Entries with lists of entries (like VEP) can be accessed in the following order:
    vcf_dict['INFO']['CSQ']['<Transcript_ID>']['<VEP_COLUMN_NAME>' (defined by format in header)]

    Args:
        vcf_file (str): path to VCF file

    Returns:
        list of dict: nested dictionary containing information from each VCF line

    Examples:
        vcf_dict = get_vcf_dict(vcf_file_lines)
        vcf_dict['REF']  # Get REF field
        vcf_dict['INFO']['MY_TAG']  # get tag from INFO field
        vcf_dict['INFO']['CSQ']['NM_0000353.3']['CONSEQUENCE']  # get nested data from tags that contain lists (like VEP)

        Note that VEP lists annotations for a number of features in the same tag. To access each I have indexed results
        by the feature IDs. If you do not know the feature ID, can iterate over all values rather than access via keys.
    """

    vcf_file_lines = [x.strip() for x in open(vcf_file, 'r')]
    return get_vcf_dict(vcf_file_lines)


def get_vcf_dict(vcf_file_lines):
    """
    Returns a list, where each each element is a nested dictionary providing access to progressively lower-level
    information from a VCF file.

    The first level is the VCF column name (all caps). The info column is a nested dictionary where tag names can be
    accessed. Entries with lists of entries (like VEP) can be accessed in the following order:
    vcf_dict['INFO']['CSQ']['<Transcript_ID>']['<VEP_COLUMN_NAME>' (defined by format in header)]

    Args:
        vcf_file_lines (list of str): list of lines from VCF file in original order

    Returns:
        list of dict: nested dictionary containing information from each VCF line

    Examples:
        vcf_dict = get_vcf_dict(vcf_file_lines)
        vcf_dict['REF']  # Get REF field
        vcf_dict['INFO']['MY_TAG']  # get tag from INFO field
        vcf_dict['INFO']['CSQ']['NM_0000353.3']['CONSEQUENCE']  # get nested data from tags that contain lists (like VEP)

        Note that VEP lists annotations for a number of features in the same tag. To access each I have indexed results
        by the feature IDs. If you do not know the feature ID, can iterate over all values rather than access via keys.

    """

    info_formats = _get_info_formats(vcf_file_lines)

    result = []
    for line in vcf_file_lines:
        line.strip()
        if not line.startswith(VCF_HEADER_PREFIX):
            columns = line.split(VCF_DELIMTER)
            column_names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
            vcf_dict = dict(zip(column_names, columns))

            vcf_dict['INFO'] = _get_info_column_dict(vcf_dict['INFO'], info_formats)
            result.append(vcf_dict)
    return result


def get_vcf_header_lines(vcf_file):
    """
    Return a list with the header lines from a VCF file.

    Args:
        vcf_file: path to VCF file

    Returns:
        list of str: list of header lines from VCF file

    """

    def stop_loop():
        raise StopIteration

    return list(x.strip("\n") if x.startswith(VCF_HEADER_PREFIX) else stop_loop() for x in open(vcf_file, 'r'))


def _get_info_column_dict(info_text, info_formats):
    """
    Returns nested dict containing information from the INFO column of VCF. See get_vcf_dict for more info. This is helper
    function to construct the dict nested at vcf_dict['INFO']

    Args:
        info_text (str): text from INFO column of VCF
        info_formats (dict): dict mapping INFO column tag IDs to respective column names

    Returns:
        Returns nested dict containing information from the INFO tag (from VEP) in info column of VCF.

    """

    info_dict = dict([x.split('=') for x in info_text.split(';')])

    for format in info_formats:
        # Further parse tags that define lists of entries like VEP
        if len(info_formats[format]) > 1:
            # Tags with list formats require different parsing
            info_dict[format] = _get_info_tag_dict(info_dict[format], info_formats[format])

    return info_dict


def _get_info_tag_dict(tag_text, info_format):
    """
    Returns a list of nested dicts given a tag from an info column of VCF. See get_vcf_dict for
    more info. This is helper function to construct the dict nested at vcf_dict['INFO'][<tag_name>] when a tag is
    in a VEP-like format.

    Args:
        tag_text (str): text from tag in INFO column of VCF that has VEP-like format
        info_formats (dict): list of column names for list of entries the tag's values

    Returns:
        Returns nested dict containing information from the CSQ tag (from VEP) in info column of VCF.

    """
    tag_values = [x.split('|') for x in tag_text.split(',')]

    result = []
    for items in tag_values:
        result.append(dict(zip(info_format, items)))

    return result


def _get_info_formats(vcf_file_lines):
    """
    Returns a dictionary of INFO column tag IDs to their respective column names from format string.

    Args:
        vcf_file_lines (str): list of lines from VCF file in original order

    Returns:
        dict: a nested dictionary of INFO column tag IDs to their respective column names from format string.

    """
    tag_pattern = re.compile('ID=(.+),')
    format_pattern = re.compile('FORMAT: (.+)"')
    format_dict = {}

    for line in vcf_file_lines:

        if not line.startswith(VCF_HEADER_PREFIX):  # only process header lines
            break

        if '##INFO' in line:
            id = _get_id_string(line)
            format = _get_format_string(line)
            format = _normalize_format_string(format)

            format_dict[id] = format.split('|')
    return format_dict


def _get_id_string(info_header_line):
    """
    Returns ID string from tag INFO entry in VCF header

    Args:
        info_header_line (str): tag INFO entry in VCF header

    Returns:
        str: ID string for tag

    """
    pattern = re.compile('id=(.+),number', re.IGNORECASE)
    match = re.search(pattern, info_header_line)
    return match.group(1)


def _get_format_string(info_header_line):
    """
    Returns format string from tag INFO entry in VCF header

    Args:
        info_header_line (str): tag INFO entry in VCF header

    Returns:
        str: format string for tag

    """
    pattern = re.compile('format: (.+)"', re.IGNORECASE)
    match = re.search(pattern, info_header_line)
    return match.group(1)


def _normalize_format_string(format_string):
    """
    Returns copy of format_string converted to all uppercase and all non-alphaneumeric characters replaced with underscore

    Args:
        info_header_line (str): format string from tag INFO entry in VCF header

    Returns:
        str: copy of format_string converted to all uppercase and all non-alphaneumeric characters replaced with underscore

    """
    format_string = format_string.upper()
    pattern = re.compile('[^A-Za-z|]')
    return re.sub(pattern, '_', format_string)


def remove_malformed_fields(data_frame):
    """
    Remove garbage entries like NM_35822.3:c.=. p.(=), r.=, etc from all entries stored in a pandas dataframe. Indended
    to pre-process table before converting to VCF format.

    Args:
        data_frame: pandas data frame

    Returns:
        data_frame with all garbage entries replaced with ''

    """

    data_frame = data_frame.replace('.+\.[\(\[\<]?(?:$|=|\?|-|0)[\)\]\>]?(?:$)', '', regex=True)
    return data_frame.replace('^-$', '', regex=True)  # entries with only dashes


def _convert_to_vcf_friendly_text(data_frame):
    """
    Replace characters that are not VCF-friendly in individual fields like ',', '|', '=', ';', etc. in pandas dataframe.
    Idended to pre-process table before converting to VCF format to ensure quality of output.

    Args:
        data_frame: pandas dataframe

    Returns:
        data_frame with ',', ';', and '|' replaced with '&'. Non-vcf-friendly characters are converted to ''.

    """
    data_frame = data_frame.replace('[=\s]', '', regex=True)
    data_frame = data_frame.replace('[,;|]', '&', regex=True)

    return data_frame


def get_vcf_info_header(data_frame, tag_id, description):
    """
    Returns INFO field header for data in dataframe. Intended for use when want to include all columns in a table as
    fields in a VCF INFO tag. This generates the header text for that INFO tag. Column titles (uppercase) are used as entry names.

    Args:
        data_frame (pandas dataframe): pandas dataframe
        tag_id (str): name for this INFO tag
        description (str): description of this INFO tag

    Returns:
        str: INFO field header for tag

    """

    format_string = '|'.join(data_frame.columns).upper()
    return '##INFO=<ID=' + tag_id + ',Number=.,TYPE=String,Description="' + description + ' Format: ' + format_string + '">'


def _map_to_genomic_coordinates(hgvs_variant, remapper):
    """
    Helper function to use with pandas.DataFrame.apply. Converts HGVS entries to VCF format.

    Args:
        hgvs_variant (str): HGVS formatted variant
        remapper (leiden.remapping.VariantRemapper): VariantRemapper (accepted as parameter because creation is expensive)

    Returns:
        pandas.DataSeries: VCF representation of variant

    """
    try:
        chrom, pos, ref, alt = remapper.hgvs_to_vcf(hgvs_variant)
        return pd.Series({'CHROM': chrom, 'POS': pos, 'ID': hgvs_variant, 'REF': ref, 'ALT': alt})
    except Exception as e:
        return pd.Series({'CHROM': '.', 'POS': '.', 'ID': '.', 'REF': '.', 'ALT': '.'})


def convert_to_vcf_format(data_frame, remapper, hgvs_column, info_tag):
    """
    Converts a pandas dataframe that contains HGVS variants along with other data about the variants to a VCF representation.
    All data fields for each entry are included in the INFO field in info_tag=column1|column2| ... format.

    Args:
        data_frame (pandas.DataFrame): pandas data_frame containing HGVS format variants and other data about said variants
        remapper (leiden.remapping.VariantRemapper): leiden.remapping.VariantRemapper object
        hgvs_column (str): name of the column containing HGVS entries
        info_tag (str): name of the tag to place in the INFO column

    Returns:
        pandas.DataFrame: new dataframe containing variants in VCF format.
            Other columns are included in INFO column as described above.

    """

    data_frame = _convert_to_vcf_friendly_text(data_frame)

    vcf_format = data_frame[hgvs_column].apply(_map_to_genomic_coordinates, args=[remapper])
    info = info_tag + '=' + data_frame.apply(lambda row: '|'.join(map(str, row)), axis=1)
    vcf_format['INFO'] = info
    vcf_format['FILTER'] = '.'
    vcf_format['QUAL'] = '.'

    vcf_format.sort(['CHROM', 'POS'])
    return vcf_format

