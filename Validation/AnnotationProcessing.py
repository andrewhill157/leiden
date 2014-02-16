import re
import glob
import traceback
import collections


def map_aa_codes(code):
    """
    Remap three letter AA codes to single letter codes. Single letter codes are returned unchanged.

    @param code: Three letter or one letter code for an amino acid or stop codon.
    @type code: string
    @return: Single letter amino acid code or * character for stop codon. Values are returned in uppercase.
    @rtype: string
    @raise: Value error if the code is not a recognized amino acid or stop codon code.
    """

    code = code.upper()

    if code == '*' or code == '':
        return code
    if code != 'X' and len(code) == 1:
        return code

    try:
        # Mapping from three letter amino acid codes to one letter amino acid codes
        one_letter_aa_codes = dict(VAL='V', ILE='I', LEU='L', GLU='E', GLN='Q', ASP='D', ASN='N', HIS='H', TRP='W', PHE='F',
                                   TYR='Y', ARG='R', LYS='K', SER='S', THR='T', MET='M', ALA='A', GLY='G', PRO='P', CYS='C',
                                   X='*', XAA='*', SCY='*')  # Last three entries are different stop-codon abbreviations

        return one_letter_aa_codes[code]

    except KeyError:
        raise ValueError('Unrecognized amino acid or stop codon code: ' + code)


# TODO would be great to make this more robust
def get_vcf_info_column(line):
    """
    Given a line of text from a tab-delimited VCF file, return the info column as a list of the included tags. Tags
    within the info column are assumed to be semi-colon separated.

    @param line: single non-header line (must be a line describing a variant) from a VCF file. Columns in the VCF text
    must be tab-separated.
    @type line: string
    @return: list of all the tags contained in the info column. List of tags is assumed to be semi-colon separated.
    @rtype: list of strings
    """

    # Lines in original file are tab separated and info is in 8th column
    info_column = line.split("\t")[7]
    return info_column.split(";")


def get_tagged_entry_value(vcf_info_entries, tag):
    """
    Returns the value of the specified tagged item from an entry from the info entry of a VCF file. Assumes that tags
    contained in vcf_info_entries take the format TAG=VALUE

    @param vcf_info_entries: a list of all the entries from the info entry of a VCF file. Note that entries
    from the info column must be contained in a list, where one tag/value pair is contained per list item.
    @type vcf_info_entries: list
    @param tag: the string used to denote the tag of interest in the info column of the VCF file
    @type tag: string
    @return: the value of the entry in the info column denoted by the specified tag. The tag and the equals sign in the
    tag notation are removed. Returns empty string if tag not found.
    @rtype: string
    """

    for entry in vcf_info_entries:
        if entry.startswith(tag):
            return remove_tag(entry)

    return ''


def remove_tag(tag_value_pair):
    """
    Removes the tag from a tag/value pair from a the info column of a VCF file.

    @param tag_value_pair: a tag value pair in the format TAG=VALUE
    @type tag_value_pair: string
    @return: value string from the tag/value pair
    @rtype: string
    @raise: ValueError if the pattern TAG=VALUE is not found
    """

    search_pattern = re.compile('.+=.+')

    if re.search(search_pattern, tag_value_pair) is not None:
        return tag_value_pair.split('=', 1)[1]
    else:
        raise ValueError('Invalid tag notation, must be in format TAG=VALUE: ' + tag_value_pair)


def remove_p_dot_notation(annotation_text):
    """
    Removes the p-dot notation from a description of protein change. Accepted formats are p.change, p.(change,
    or p.[change], where change is returned and is the description of the protein change.

    @param annotation_text: p-dot notation describing the protein change
    @type annotation_text: string
    @return: annotation_text with the p-dot notation removed
    @rtype: string
    @raise: ValueError if the p-dot notation is not in one of the expected formats.
    """

    search_pattern = re.compile('[pP]\.[\(\[]?([^\)\]]+)[\)\]]?')
    match = re.search(search_pattern, annotation_text)

    if match is not None:
        return match.group(1)
    elif annotation_text == '-':
        return annotation_text
    else:
        raise ValueError('Invalid p-dot notation. Should be p.change, p.(change) or p.[change]: ' + annotation_text)


def get_laa_change(vcf_info_column_list):
    """
    Returns the before and after amino acid or stop codon values from the HGVS notation protein change entry from LOVD,
    denoted by the tag LAA_CHANGE in the INFO column of a VCF file.

    @param vcf_info_column_list: list of tag/value pairs from entry in the info column of VCF file. Each tag must take
    the format TAG=VALUE.
    @type vcf_info_column_list: list
    @return: before and after amino acids or stop codons from the HGVS notation protein change entry from LOVD.
    @rtype: named tuple containing string elements named before and after
    @raise: ValueError if annotation info does not match the expected format.
    """

    laa_change = collections.namedtuple('laa_change', 'before after')

    laa_change_tag = "LAA_CHANGE"
    laa_change_value = get_tagged_entry_value(vcf_info_column_list, laa_change_tag)
    laa_change_value = remove_p_dot_notation(laa_change_value)

    # TODO - temporary fix. Ignore indels that have been accidentally included
    if 'ins' in laa_change_value or 'del' in laa_change_value or 'fs' in laa_change_value:
        raise ValueError('Indel Detected: ' + laa_change_value)


    # Protein change values that indicate that there is no protein change
    denotes_no_aa_change = ['-', '=', '?']

    if laa_change_value in denotes_no_aa_change:
        return laa_change('', '')
    else:
        pattern = re.compile('([A-Za-z\*\?]{1,3})[\d+]+([A-Za-z\*]{1,3})[A-Za-z\*\?]*\d*')
        match = re.search(pattern, laa_change_value)

        if match is not None:
            return laa_change(match.group(1), match.group(2))
        else:
            raise ValueError('Unexpected laa_change format: ' + laa_change_value)


def get_aa_change(vcf_info_column_list):
    """
    Return the before and after amino acid or stop codon change from VEP, as denoted by the AA_CHANGE tag in the INFO
    column of a VCF file.

    @param vcf_info_column_list: list of tag/value pairs from entry in the info column of VCF file. Each tag must take
    the format TAG=VALUE.
    @type vcf_info_column_list: list
    @return: before and after amino acids or stop codons from the VEP AA_CHANGE change entry. Before and after will be
    empty strings if no match was found, indicating no predicted change.
    @rtype: named tuple containing string elements named before and after
    """

    aa_change = collections.namedtuple('aa_change', 'before after')

    aa_change_tag = "AA_CHANGE"
    raw_annotation = get_tagged_entry_value(vcf_info_column_list, aa_change_tag)
    raw_annotation = raw_annotation.split(",")

    for change in raw_annotation:
        search_pattern = re.compile('([A-Za-z\*]+)/([A-Za-z\*]+)?')
        match = re.search(search_pattern, change)

        if match is not None:
            return aa_change(match.group(1), match.group(2))
        else:
            # Check for synonymous changes which are a single letter
            search_pattern = '([A-Za-z\*])'
            match = re.search(search_pattern, change)

            if match is not None:
                # Return with before and after being the same amino acid
                return aa_change(match.group(1), match.group(1))

    return aa_change('', '')  # no entry found


def is_concordant_annotation(laa_change, aa_change):
    """
    Given the laa_change and aa_change tagged entries (the Protein change entry from LOVD and the predicted changes from
    VEP), determine whether the two annotations are concordant.

    @param laa_change: before and after amino acids or stop codons from the HGVS notation protein change entry from LOVD.
    @type laa_change: named tuple containing string elements named before and after
    @param aa_change: before and after amino acids or stop codons from the VEP AA_CHANGE change entry.
    @type aa_change: named tuple containing string elements named before and after
    @return: true if the laa_change and aa_change are concordant, false otherwise
    @rtype: bool
    """
    laa_change_before = map_aa_codes(laa_change.before)
    laa_change_after = map_aa_codes(laa_change.after)
    aa_change_before = map_aa_codes(aa_change.before)
    aa_change_after = map_aa_codes(aa_change.after)

    matching_annotation = laa_change_before == aa_change_before and laa_change_after == aa_change_after
    synonymous_mutation = laa_change_before == laa_change_after and aa_change_before == aa_change_after

    return matching_annotation or synonymous_mutation


def get_severe_impact(vcf_info_column_list):
    """
    Return the value of the SEVERE_IMPACT entry from the INFO column of the VCF file.

    @param vcf_info_column_list: list of tag/value pairs from entry in the info column of VCF file. Each tag must take
    the format TAG=VALUE.
    @type vcf_info_column_list: list
    @return: value of the SEVERE_IMPACT entry from the INFO column of the VCF file.
    @rtype: string
    """
    severe_impact_flag = "SEVERE_IMPACT"
    return get_tagged_entry_value(vcf_info_column_list, severe_impact_flag)


def get_overall_26K_allele_frequency(vcf_info_column_list):
    """
    Given a list of INFO column entries for a mutation from VCF file, return the overall allele frequency within the
    26K data.

    @param vcf_info_column_list: list of tag/value pairs from entry in the info column of VCF file. Each tag must take
    the format TAG=VALUE.
    @type vcf_info_column_list: list
    @return: overall allele frequency within the 26K data. Returns 0 if not found within 26K data.
    @rtype: float
    """

    ac_flag = 'AC_MAC26K'
    an_flag = 'AN_MAC26K'

    ac = get_tagged_entry_value(vcf_info_column_list, ac_flag)
    an = get_tagged_entry_value(vcf_info_column_list, an_flag)

    if ac == '' or an == '':
        return 0
    else:
        ac = float(ac)
        an = float(an)
        return (ac/an)*100

