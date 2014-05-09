import re


def is_concordant(variant):
    try:
        lovd_protein_change = variant['INFO']['LOVD']['PROTEIN_CHANGE'].lower()
        vep_transcripts = variant['INFO']['CSQ']
    except KeyError:
        return False

    for transcript_id in vep_transcripts:

        lovd_protein_change = re.sub('xaa', '*', lovd_protein_change)
        lovd_protein_change = re.sub('x', '*', lovd_protein_change)
        lovd_protein_change = remove_p_dot_notation(lovd_protein_change)

        vep_hgvs_protein_change = vep_transcripts[transcript_id]['HGVSP']

        if lovd_protein_change and vep_hgvs_protein_change:
            vep_hgvs_protein_change = vep_hgvs_protein_change.split(':')[1].lower()
            vep_hgvs_protein_change = re.sub('ter', '*', vep_hgvs_protein_change)
            vep_hgvs_protein_change = remove_p_dot_notation(vep_hgvs_protein_change)

            if lovd_protein_change == vep_hgvs_protein_change:
                return True

    return False


def get_ucsc_location_link(chromosome_number, start_coordinate, end_coordinate):
    """
    Returns link to relevant range in the UCSC genome browser. All parameters must be in valid range.

    Args:
        chromosome_number (str): the chromosome number of the range to link to
        start_coordinate (str): the start coordinate of range to link to
        end_coordinate (str): the end coordinate of range to link to

    Returns:
        str: URL link to display region in UCSC genome browser

    """
    return ''.join(['http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr', chromosome_number, '%3A',
                    start_coordinate, '-', end_coordinate])


def remove_p_dot_notation(annotation_text):
    """
    Removes the p-dot notation from a description of protein change. Accepted formats are p.change, p.(change,
    or p.[change], where change is returned and is the description of the protein change.

    Args:
        annotation_text (str): p-dot notation describing the protein change

    Returns:
        str: Annotation_text with the p-dot notation removed

    Exceptiom:
        ValueError: if the p-dot notation is not in one of the expected formats.

    """

    search_pattern = re.compile('[p]\.[\(\[]?([^\)\]]+)[\)\]]?', re.IGNORECASE)
    match = re.search(search_pattern, annotation_text)

    if match:
        return match.group(1)
    elif annotation_text == '-':
        return annotation_text
    else:
        return annotation_text