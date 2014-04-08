"""
Andrew Hill
MacArthur Lab - 2014

Outputs error logs and validation statistics on VEP and 26K annotated VCF files. Meant to validate variants converted
from HGVS to VCF format.

For help, execute: python process_annotated_mutations.py --help
"""

import argparse
import os
from macarthur_core.io import file_io
import re
from matplotlib import pyplot

CHROMOSOME_NUMBER_INDEX = 0
COORDINATE_INDEX = 1
REF_INDEX = 3
ALT_INDEX = 4
INFO_COLUMN_INDEX = 7


def has_vep_aa_change(variant_vcf_row):
    try:
        info_column = variant_vcf_row[INFO_COLUMN_INDEX]
        get_unique_tagged_entry_values(info_column, 'AA_CHANGE')
        return True
    except ValueError:
        return False


def has_lovd_aa_change(variant_vcf_row):
    try:
        info_column = variant_vcf_row[INFO_COLUMN_INDEX]
        laa_change_value = get_unique_tagged_entry_values(info_column, 'LAA_CHANGE')[0]
        laa_change_value = remove_p_dot_notation(laa_change_value)
    except:
        return False

    # Protein change values that indicate that there is no reported protein change
    denotes_no_aa_change = ['-', '=', '?', '0']

    return laa_change_value not in denotes_no_aa_change


def get_severe_impact(variant_vcf_row):
    info_column = variant_vcf_row[INFO_COLUMN_INDEX]

    return get_unique_tagged_entry_values(info_column, 'SEVERE_IMPACT')[0]


def is_concordant(variant_vcf_row):

    lovd_before, lovd_after = get_laa_change(variant_vcf_row)
    lovd_before = map_aa_codes(lovd_before)
    lovd_after = map_aa_codes(lovd_after)

    vep_aa_change = get_vep_aa_change(variant_vcf_row)

    for entry in vep_aa_change:
        vep_before, vep_after = entry

        if lovd_before == vep_before and lovd_after == vep_after:
            return True

    return False


def is_concordant_splice_mutation(variant_vcf_row):
    info_column = variant_vcf_row[INFO_COLUMN_INDEX]
    hgvs = get_unique_tagged_entry_values(info_column, 'HGVS')[0]
    ref = variant_vcf_row[REF_INDEX]

    pattern = re.compile('([+-]\d)([A-Z])>[A-Z]', re.IGNORECASE)
    match = re.search(pattern, hgvs)

    if match:
        lovd_splice_position = match.group(1)
        lovd_splice_base = match.group(2)
    else:
        raise ValueError('Unexpected splice variant pattern: ' + hgvs)

    vep_splice_position_list = get_unique_tagged_entry_values(info_column, 'SPLICE_POS')
    conserved_splice_combinations = ['+1G', '+2T', '-2A', '-1G']

    for vep_position in vep_splice_position_list:
        if not vep_position.startswith('-'):
            vep_position = '+' + vep_position

        if vep_position + ref == lovd_splice_position + lovd_splice_base and lovd_splice_position + lovd_splice_base in conserved_splice_combinations:
            return True
    return False


def is_concordant_frameshift_mutation(variant_vcf_row):
    return False


def is_concordant_inframe_codon_loss(variant_vcf_row):
    return False


def plot_allele_frequency_histogram(frequency_list):
    """
    Generate histogram plots of allele frequency data
    """

    frequency_list = [x if x < 10 else 10 for x in frequency_list]
    p = pyplot.hist(frequency_list, 100, facecolor='blue', alpha=0.75)
    vertical_line = pyplot.axvline(x=0.5, color='red', ls='dashed')
    
    
def get_ucsc_location_link(chromosome_number, start_coordinate, end_coordinate):
    """
    Returns link to relevant range in the UCSC genome browser. All parameters must be in valid range.

    @param chromosome_number: the chromosome number of the range to link to
    @type chromosome_number: string
    @param start_coordinate: the start coordinate of range to link to
    @type start_coordinate: string
    @param end_coordinate: the end coordinate of range to link to
    @type end_coordinate: string
    @return: URL link to display region in UCSC genome browser
    @rtype: string
    """
    return ''.join(['http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr', chromosome_number, '%3A',
                    start_coordinate, '-', end_coordinate])


def map_aa_codes(code):
    """
    Remap three letter AA codes to single letter codes. Single letter codes are returned unchanged.

    @param code: Three letter or one letter code for an amino acid or stop codon.
    @type code: string
    @return: Single letter amino acid code or * character for stop codon. Values are returned in uppercase.
    @rtype: string
    @raise: Value error if the code is not a recognized three-letter amino acid or stop codon code.
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
                                   X='*', XAA='*', SCY='*',  # stop codon abbreviations
                                   DEL='DEL')  # Deletion

        return one_letter_aa_codes[code]

    except KeyError:
        raise ValueError('Unrecognized amino acid or stop codon code: ' + code)


def get_unique_tagged_entry_values(info_column, tag, tagged_entry_delimiter=';', value_delimiter=','):
    """
    Returns the unique values of the specified tagged item from an entry from the info column of a VCF file.
    Assumes that tags contained in vcf_info_entries take the format TAG=VALUE, where value is assumed to be a list of
    values.

    @param info_column: text from the info entry of a VCF file. Note that entries from the info column must be contained
    in a list, where tag/value pairs are used <tag>=<value>;<tag1>=<value1>.
    @type info_column: string
    @param tag: the string used to denote the tag of interest in the info column of the VCF file
    @type tag: string
    @param tagged_entry_delimiter: delimiter for sequence of tag-value pairs
    @type tagged_entry_delimiter: string
    @param value_delimiter: delimiter for sequence of values associated with tagged entries
    @type value_delimiter: string
    @return: list of unique values for the entry in the info column denoted by the specified tag. All empty entries
    (empty strings or - characters) will not be included in results. Therefore, an empty list is returned if tag is
    present but no non-empty value is found.
    @rtype: list
    @raises: ValueError if tag not found in info_column
    """

    pattern = re.compile('(?:^|' + tagged_entry_delimiter + ')' + tag + '=([^;]*)(?:' + tagged_entry_delimiter + '|$)', re.IGNORECASE)
    match = re.search(pattern, info_column)

    if match:
        # Make a list of all non-empty entries
        result = [value for value in match.group(1).split(value_delimiter) if value != '-' and value != '']
        return list(set(result))
    else:
        raise ValueError('The requested tag or tag entry delimiter was not found: ' + tag)


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

    search_pattern = re.compile('[p]\.[\(\[]?([^\)\]]+)[\)\]]?', re.IGNORECASE)
    match = re.search(search_pattern, annotation_text)

    if match:
        return match.group(1)
    elif annotation_text == '-':
        return annotation_text
    else:
        raise ValueError('Invalid p-dot notation. Should be p.change, p.(change) or p.[change]: ' + annotation_text)


def get_laa_change(variant_vcf_row):
    """
    Returns the before and after amino acid or stop codon values from the HGVS notation protein change entry from lovd,
    denoted by the tag LAA_CHANGE in the INFO column of a VCF file.

    @param vcf_info_column_list: list of tag/value pairs from entry in the info column of VCF file. Each tag must take
    the format TAG=VALUE.
    @type vcf_info_column_list: list
    @return: before and after amino acids or stop codons from the HGVS notation protein change entry from lovd.
    @rtype: named tuple containing string elements named before and after
    @raise: ValueError if annotation info does not match the expected format.
    """
    info_column = variant_vcf_row[INFO_COLUMN_INDEX]

    laa_change_tag = 'LAA_CHANGE'

    laa_change_value = get_unique_tagged_entry_values(info_column, laa_change_tag)[0]

    laa_change_value = remove_p_dot_notation(laa_change_value)

    pattern = re.compile('([A-Za-z\*\?]{1,3})[\d]+([A-Za-z\*]{1,3})[A-Za-z\*\?]*\d*')
    match = re.search(pattern, laa_change_value)

    if match is not None:
        return match.group(1), match.group(2)
    else:
        raise ValueError('Unexpected laa_change format: ' + laa_change_value)


def get_vep_aa_change(variant_vcf_row):
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
    info_column = variant_vcf_row[INFO_COLUMN_INDEX]
    aa_change_tag = 'AA_CHANGE'
    raw_annotation = get_unique_tagged_entry_values(info_column, aa_change_tag)
    aa_change_list = []

    for change in raw_annotation:
        search_pattern = re.compile('([A-Za-z\*]+)/([A-Za-z\*]+)?')
        match = re.search(search_pattern, change)

        if match is not None:
            aa_change_list.append((match.group(1), match.group(2)))
        else:
            # Check for synonymous changes which are a single letter
            search_pattern = '([A-Za-z\*])'
            match = re.search(search_pattern, change)

            if match:
                # Return with before and after being the same amino acid
                aa_change_list.append((match.group(1), match.group(1)))

    return aa_change_list


def get_26K_allele_frequency(population_set, variant_vcf_row):
    """
    Given a list of INFO column entries for a mutation from VCF file, return the overall allele frequency within the
    26K data.

    @param population_set: tag for population set within 26K data set. MAC26K is the overall set, SAS is south asian, etc.
    These should match AC_<population_set> and AN_<population_set> tags in vcf_info_column_list.
    @type population_set: string
    @param vcf_info_column_list: list of tag/value pairs from entry in the info column of VCF file. Each tag must take
    the format TAG=VALUE.
    @type vcf_info_column_list: list
    @return: allele frequency within specific population set in 26K data. Returns 0 if not found within 26K data.
    @rtype: float
    """

    info_column = variant_vcf_row[INFO_COLUMN_INDEX]

    ac_flag = 'AC_' + population_set
    an_flag = 'AN_' + population_set

    try:
        ac = get_unique_tagged_entry_values(info_column, ac_flag)[0]
        an = get_unique_tagged_entry_values(info_column, an_flag)[0]
    except:
        return 0

    ac = float(ac)
    an = float(an)
    return (ac/an)*100

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Output validation statitics and error logs when processing VEP and AC '
                                                 'annotated VCF files. Outputs processing_errors.log and '
                                                 'discordant_annotations.log, which contain information on variants that '
                                                 'could not be validated due to syntax errors and variants that not '
                                                 'correctly validate respectively.')
    
    group = parser.add_argument_group()
    group.add_argument('-f', '--file_names', required=True, help='File containing full paths to the VCF files to be processed.')
    group.add_argument('-o', '--output_directory', default='.', help='Output directory for .log files.')
    group.add_argument('-p', '--plots', action='store_true', help='Set to enable generation of summary plots.')
    args = parser.parse_args()
    
    # Make the output directory if does not already exist
    output_directory = args.output_directory
    
    if not os.path.exists(output_directory):
        os.mkdir(args.output_directory)
    
    
    # Counts of various categories of variants
    no_predicted_aa_change_count = 0
    missing_lovd_aa_change_count = 0
    concordant_annotation_count = 0
    discordant_annotation_count = 0
    error_count = 0
    not_processed_count = 0
    
    hgmd_site_count = 0
    hgmd_mutation_count = 0
    
    overall_high_26K_frequency_count = []
    eur_high_26K_frequency_count = []
    amr_high_26K_frequency_count = []
    afr_high_26K_frequency_count = []
    sas_high_26K_frequency_count = []
    eas_high_26K_frequency_count = []
    
    overlap_26K_count = 0
    
    total_mutation_count = 0
    
    no_vep_aa_change = []
    processing_errors = []
    discordant_annotations = []
    
    with open(args.file_names, 'r') as file_list:
        files_to_process = file_list.read().splitlines()
    
    for file in files_to_process:
        # Extract all non-header lines from VCF file
        vcf_file = file_io.read_table_from_file(file)
        variants = [x for x in vcf_file if not x[0].startswith('#')]
    
        for variant_row in variants:

            total_mutation_count += 1

            # Extract info from line of VCF
            chromosome_number = variant_row[CHROMOSOME_NUMBER_INDEX]
            coordinate = variant_row[COORDINATE_INDEX]
            ref = variant_row[REF_INDEX]
            alt = variant_row[ALT_INDEX]
            info = variant_row[INFO_COLUMN_INDEX]
    
            # Get link to relevant location in UCSC genome browser
            viewing_interval = 25
            if coordinate != '.':
                ucsc_link = get_ucsc_location_link(chromosome_number,
                                                   str(int(coordinate) - viewing_interval),
                                                   str(int(coordinate) + viewing_interval))
    
            # Check allele frequency
            overall_allele_frequency = get_26K_allele_frequency('MAC26K', variant_row)
            eur_allele_frequency = get_26K_allele_frequency('EUR', variant_row)
            amr_allele_frequency = get_26K_allele_frequency('AMR', variant_row)
            afr_allele_frequency = get_26K_allele_frequency('AFR', variant_row)
            sas_allele_frequency = get_26K_allele_frequency('SAS', variant_row)
            eas_allele_frequency = get_26K_allele_frequency('EAS', variant_row)
    
            HIGH_FREQUENCY_THRESHOLD = 0
            if overall_allele_frequency > HIGH_FREQUENCY_THRESHOLD:
                overall_high_26K_frequency_count.append(overall_allele_frequency)
            if eur_allele_frequency > HIGH_FREQUENCY_THRESHOLD:
                eur_high_26K_frequency_count.append(eur_allele_frequency)
            if amr_allele_frequency > HIGH_FREQUENCY_THRESHOLD:
                amr_high_26K_frequency_count.append(amr_allele_frequency)
            if afr_allele_frequency > HIGH_FREQUENCY_THRESHOLD:
                afr_high_26K_frequency_count.append(afr_allele_frequency)
            if sas_allele_frequency > HIGH_FREQUENCY_THRESHOLD:
                sas_high_26K_frequency_count.append(sas_allele_frequency)
            if eas_allele_frequency > HIGH_FREQUENCY_THRESHOLD:
                eas_high_26K_frequency_count.append(eas_allele_frequency)
    
            if overall_allele_frequency > 0:
                overlap_26K_count += 1
    
            # Check HGMD Overlap
            try:
                HGMD_SITE_TAG = 'HGMD_SITE'
                hgmd_site = get_unique_tagged_entry_values(info, 'HGMD_SITE')
                hgmd_site_count += 1
            except ValueError:
                pass

            try:
                HGMD_MUTATION_TAG = 'HGMD_MUT'
                hgmd_mutation = get_unique_tagged_entry_values(info, 'HGMD_MUT')
                hgmd_mutation_count += 1
            except ValueError:
                pass
    
            # Check concordance
            if has_vep_aa_change(variant_row) and not has_lovd_aa_change(variant_row):
                # If VEP provided AA change, LOVD should have reported AA change. Cannot validate.
                missing_lovd_aa_change_count += 1
            else:
                try:
                    severe_impact = get_severe_impact(variant_row)

                    if has_vep_aa_change(variant_row):
                        validation_function = is_concordant
                    elif severe_impact == 'INFRAME_CODON_LOSS':
                        not_processed_count += 1
                        validation_function = is_concordant_inframe_codon_loss
                    elif severe_impact == 'SPLICE_ACCEPTOR_VARIANT' or severe_impact == 'SPLICE_DONOR_VARIANT':
                        validation_function = is_concordant_splice_mutation
                    elif severe_impact == 'FRAMESHIFT_VARIANT':
                        not_processed_count += 1
                        validation_function = is_concordant_frameshift_mutation
                    else:
                        # Variant is not a category that we can process right now
                        not_processed_count += 1
                        no_vep_aa_change.append([severe_impact] + variant_row)
                        validation_function = is_concordant_frameshift_mutation

                    concordant = validation_function(variant_row)

                    if concordant:
                        concordant_annotation_count += 1
                    else:
                        discordant_annotation_count += 1
                        if has_vep_aa_change(variant_row):
                            discordant_annotations.append([severe_impact, get_unique_tagged_entry_values(info, 'AA_CHANGE')[0], get_unique_tagged_entry_values(info, 'LAA_CHANGE')[0]])
                        else:
                            discordant_annotations.append([severe_impact, get_unique_tagged_entry_values(info, 'LAA_CHANGE')[0]])
                except Exception as e:
                    error_count += 1
                    processing_errors.append([str(e)])
    
    # Write out entries with no VEP AA change for further examination
    no_vep_aa_change_file_name = os.path.join(output_directory, 'no_vep_aa_change.log')
    file_io.write_table_to_file(no_vep_aa_change_file_name, no_vep_aa_change)
    
    # Write out errors and discordant mutations to file
    processing_errors.insert(0, ['file',
                                 'error',
                                 'hgvs',
                                 'protein',
                                 'severe_impact'])
    
    processing_errors_file_name = os.path.join(output_directory, 'processing_errors.log')
    file_io.write_table_to_file(processing_errors_file_name, processing_errors)
    
    discordant_annotations.insert(0, ['file',
                                      'hgvs',
                                      'chromosome',
                                      'coordinate',
                                      'ref',
                                      'alt',
                                      'ucsc',
                                      'protein',
                                      'aa_change_lovd',
                                      'aa_change_vep',
                                      'severe_impact'])
    discordant_annotations_file_name = os.path.join(output_directory, 'discordant_annotations.log')
    file_io.write_table_to_file(discordant_annotations_file_name, discordant_annotations)
    
    # Print results
    print 'Total Mutations: ' + str(total_mutation_count)
    print 'Missing LOVD AA Change (Unable to Validate): ' + str(missing_lovd_aa_change_count)
    print 'Not Processed: ' + str(not_processed_count)
    print 'Concordant Annotations: ' + str(concordant_annotation_count)
    print 'Discordant Annotations: ' + str(discordant_annotation_count)
    print 'Errors: ' + str(error_count)
    print ''
    print 'HGMD Sites: ' + str(hgmd_site_count)
    print 'HGMD Mutations: ' + str(hgmd_mutation_count)

    if args.plots:
        # Plot allele frequency histograms
        common_axes = pyplot.subplot(111)
        common_axes.set_xlabel('Allele Frequency (%)')
        common_axes.set_ylabel('Counts')

        pyplot.subplot(231)
        pyplot.title('Overall')
        plot_allele_frequency_histogram(overall_high_26K_frequency_count)

        pyplot.subplot(232)
        pyplot.title('European')
        plot_allele_frequency_histogram(eur_high_26K_frequency_count)

        pyplot.subplot(233)
        pyplot.title('American')
        plot_allele_frequency_histogram(amr_high_26K_frequency_count)

        pyplot.subplot(234)
        pyplot.title('African')
        plot_allele_frequency_histogram(afr_high_26K_frequency_count)

        pyplot.subplot(235)
        pyplot.title('East Asian')
        plot_allele_frequency_histogram(eas_high_26K_frequency_count)

        pyplot.subplot(236)
        pyplot.title('South Asian')
        plot_allele_frequency_histogram(sas_high_26K_frequency_count)

        pyplot.tight_layout(h_pad=2.0)
        pyplot.show()
