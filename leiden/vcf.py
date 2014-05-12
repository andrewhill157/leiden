import re
import pandas as pd
from _ordereddict import ordereddict

VCF_HEADER_PREFIX = '#'
VCF_DELIMITER = '\t'
FORMAT_DELIMITER = '|'


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
    Replace characters that are not VCF-friendly in individual fields like ',', FORMAT_DELIMITER, '=', ';', etc. in pandas dataframe.
    Idended to pre-process table before converting to VCF format to ensure quality of output.

    Args:
        data_frame: pandas dataframe

    Returns:
        data_frame with ',', ';', and FORMAT_DELIMITER replaced with '&'. Non-vcf-friendly characters are converted to ''.

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

    format_string = FORMAT_DELIMITER.join(data_frame.columns).upper()
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
    info = info_tag + '=' + data_frame.apply(lambda row: FORMAT_DELIMITER.join(map(str, row)), axis=1)
    vcf_format['INFO'] = info
    vcf_format['FILTER'] = '.'
    vcf_format['QUAL'] = '.'

    vcf_format.sort(['CHROM', 'POS'])
    return vcf_format


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

    return list(x.strip("\n") if x.startswith(VCF_HEADER_PREFIX) else stop_loop() for x in vcf_file)


class VCFReader():
    def __init__(self, file_object):
        self.file_object = file_object
        self.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        self.header_lines = get_vcf_header_lines(self.file_object)
        self.vcf_format, self.infos = self.parse_vcf_header(self.header_lines)

    def parse_vcf_header(self, vcf_header_lines):
        infos = {}
        for line in vcf_header_lines:
            if '##fileformat' in line:
                file_format = line
            elif '##INFO' in line:
                id = re.search('ID=([^,]+)\,', line).group(1)
                number = re.search('Number=([^,]+),', line).group(1)
                data_type = re.search('Type=([^,]+),', line).group(1)
                description = re.search('Description="(.+)"', line).group(1)
                infos[id] = {'number': number, 'type': data_type, 'description': description}

                if 'Format' in description:
                    tag_format = re.search('Format: (.+)', description).group(1)
                    tag_format = _normalize_format_string(tag_format)
                    infos[id]['format'] = tag_format.split(FORMAT_DELIMITER)

        return file_format, infos

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        line = self.file_object.next()
        line = line.strip()
        columns = line.split(VCF_DELIMITER)
        vcf_dict = ordereddict(zip(self.columns, columns))

        vcf_dict['INFO'] = self.get_info_dict(vcf_dict['INFO'])
        return VCFLine(vcf_dict)

    def get_info_dict(self, info_text):
        tags = [x.split('=') for x in info_text.split(';')]
        info_dict = ordereddict(tags)

        for tag in info_dict:
            if 'format' in self.infos[tag]:
                result = []
                for entry in info_dict[tag].split(','):
                    result.append(ordereddict(zip(self.infos[tag]['format'], entry.split(FORMAT_DELIMITER))))
                info_dict[tag] = result
        return info_dict


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


class VCFLine:

    def __init__(self, vcf_line):
        self.variant = vcf_line

    def __getitem__(self, var, **args):
        return self.variant[var]

    def __str__(self):
        info_column = self.variant['INFO']

        values = info_column.values()

        for i,tag_value in enumerate(values):

            if isinstance(tag_value, list):

                inner_entries = []
                for item in tag_value:
                    inner_entries.append(FORMAT_DELIMITER.join(item.values()))

                values[i] = ','.join(inner_entries)

        keys = info_column.keys()
        info = ';'.join(['%s=%s' % (keys[i], value) for i,value in enumerate(values)])
        output_string = self.variant.values()
        output_string[-1] = info
        return '\t'.join(output_string)

