#######################################################################################################################
# Tests for has_vep_aa_change
#######################################################################################################################
from leiden.bin import validate_annotated_vcfs


def test_has_vep_aa_change_with_standard_input():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'AA_CHANGE=p.(Tyr657Gly)'
    result = True
    assert_equals(validate_annotated_vcfs.has_vep_aa_change(input), result)


def test_has_vep_aa_change_with_no_vep_aa_change():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'LAA_CHANGE=p.(Tyr657Gly)'
    result = False
    assert_equals(validate_annotated_vcfs.has_vep_aa_change(input), result)


#######################################################################################################################
# Tests for has_lovd_aa_change
#######################################################################################################################
def test_has_lovd_aa_change_with_standard_input():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'LAA_CHANGE=p.(Tyr657Gly)'
    result = True
    assert_equals(validate_annotated_vcfs.has_lovd_aa_change(input), result)


def test_has_vep_aa_change_with_no_aa_change():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'AA_CHANGE=p.(Tyr657Gly)'
    result = False
    assert_equals(validate_annotated_vcfs.has_lovd_aa_change(input), result)


def test_has_vep_aa_change_with_invalid_aa_change():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'AA_CHANGE=p.(=)'
    result = False
    assert_equals(validate_annotated_vcfs.has_lovd_aa_change(input), result)


#######################################################################################################################
# Tests for get_severe_impact
#######################################################################################################################
def test_get_severe_impact_with_standard_input():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'SEVERE_IMPACT=NONSENSE_VARIANT'
    result = 'NONSENSE_VARIANT'
    assert_equals(validate_annotated_vcfs.get_severe_impact(input), result)


def test_get_severe_impact_with_multiple_categories():
    # Even if there are multiple categories, only the first one should be returned
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'SEVERE_IMPACT=NONSENSE_VARIANT,OTHER_CATEGORY'
    result = 'NONSENSE_VARIANT'
    assert_equals(validate_annotated_vcfs.get_severe_impact(input), result)


def test_get_severe_impact_with_multiple_categories():
    # Even if there are multiple categories, only the first one should be returned
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'NO_SEVERE_IMPACT_TAG=NONE'
    assert_raises(ValueError, validate_annotated_vcfs.get_severe_impact, input)


#######################################################################################################################
# Tests for is_concordant
#######################################################################################################################
def test_is_concordant_with_standard_input():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'LAA_CHANGE=p.(Gly456Tyr);AA_CHANGE=G/Y'
    result = True
    assert_equals(validate_annotated_vcfs.is_concordant(input), result)


def test_is_concordant_with_discordant_input():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'LAA_CHANGE=p.(Gly456Tyr);AA_CHANGE=G/W'
    result = False
    assert_equals(validate_annotated_vcfs.is_concordant(input), result)


def test_is_concordant_with_one_correct_option_in_list():
    # If there are multiple annotations (different transcripts), return true if any are concordant
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'LAA_CHANGE=p.(Gly456Tyr);AA_CHANGE=G/W,G/Y'
    result = True
    assert_equals(validate_annotated_vcfs.is_concordant(input), result)


#######################################################################################################################
# Tests for is_concordant_splice_mutation
#######################################################################################################################

# First four tests check for conserved +1G,+2T,-1G,-2A patterns at donor and acceptor sites
def test_is_concordant_splice_mutation_with_plus_1_g():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'HGVS=NM_00012.2:c.990+1G>T;SPLICE_POS=1'
    input[validate_annotated_vcfs.REF_INDEX] = 'G'
    result = True
    assert_equals(validate_annotated_vcfs.is_concordant_splice_mutation(input), result)


def test_is_concordant_splice_mutation_with_plus_2_t():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'HGVS=NM_00012.2:c.990+2T>G;SPLICE_POS=2'
    input[validate_annotated_vcfs.REF_INDEX] = 'T'
    result = True
    assert_equals(validate_annotated_vcfs.is_concordant_splice_mutation(input), result)


def test_is_concordant_splice_mutation_with_minus_1_g():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'HGVS=NM_00012.2:c.990-1G>T;SPLICE_POS=-1'
    input[validate_annotated_vcfs.REF_INDEX] = 'G'
    result = True
    assert_equals(validate_annotated_vcfs.is_concordant_splice_mutation(input), result)


def test_is_concordant_splice_mutation_with_minus_2_a():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'HGVS=NM_00012.2:c.990-2A>G;SPLICE_POS=-2'
    input[validate_annotated_vcfs.REF_INDEX] = 'A'
    result = True
    assert_equals(validate_annotated_vcfs.is_concordant_splice_mutation(input), result)


# Other tests
def test_is_concordant_splice_mutation_with_discordant_site():
    # This position does not match a conserved pattern
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'HGVS=NM_00012.2:c.990-2A>G;SPLICE_POS=-2'
    input[validate_annotated_vcfs.REF_INDEX] = 'G'
    result = False
    assert_equals(validate_annotated_vcfs.is_concordant_splice_mutation(input), result)


def test_is_concordant_splice_mutation_with_discordant_site():
    # This position does not match a conserved pattern
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'HGVS=NM_00012.2:c.990-2A>G;SPLICE_POS=2'
    input[validate_annotated_vcfs.REF_INDEX] = 'A'
    result = False
    assert_equals(validate_annotated_vcfs.is_concordant_splice_mutation(input), result)

def test_is_concordant_splice_mutation_with_non_conserved_acceptor_site():
    # Other splice sites are not necessarily conserved, cannot validate
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'HGVS=NM_00012.2:c.990-100A>G;SPLICE_POS=-100'
    input[validate_annotated_vcfs.REF_INDEX] = 'A'
    result = False
    assert_equals(validate_annotated_vcfs.is_concordant_splice_mutation(input), result)


def test_is_concordant_splice_mutation_with_non_conserved_acceptor_site():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'HGVS=NM_00012.2:c.990-2A>G;SPLICE_POS=-2'
    input[validate_annotated_vcfs.REF_INDEX] = 'A'
    result = True
    assert_equals(validate_annotated_vcfs.is_concordant_splice_mutation(input), result)


#######################################################################################################################
# Tests for is_concordant_inframe_codon_loss
#######################################################################################################################

# TODO - not implemented

#######################################################################################################################
# Tests for plot_allele_frequency_histogram
#######################################################################################################################

# TODO - is not easily testable with unit tests

#######################################################################################################################
# Tests for get_ucsc_location_link
#######################################################################################################################
def test_get_ucsc_location_link_with_standard_interval():
    # Basic test to ensure URL is constructed correctly
    chromosome = '1'
    start_coordinate = '15613'
    end_coordinate = '52552'
    result = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1%3A15613-52552'

    assert_equals(validate_annotated_vcfs.get_ucsc_location_link(chromosome, start_coordinate, end_coordinate), result)


#######################################################################################################################
# Tests for map_aa_codes
#######################################################################################################################
def test_map_aa_codes_with_three_letter_code():
    # Remap three letter code to one letter code
    input = 'GLU'
    result = 'E'
    assert_equals(validate_annotated_vcfs.map_aa_codes(input), result)


def test_map_aa_code_with_lower_case():
    # Case should not matter
    input = 'met'
    result = 'M'
    assert_equals(validate_annotated_vcfs.map_aa_codes(input), result)


def test_map_aa_code_with_one_letter_code():
    # One letter codes returned unchanged
    input = 'E'
    result = 'E'
    assert_equals(validate_annotated_vcfs.map_aa_codes(input), result)


def test_map_aa_code_with_lower_case_one_letter_code():
    # One letter codes are returned as capitals regardless of input case
    input = 'k'
    result = 'K'
    assert_equals(validate_annotated_vcfs.map_aa_codes(input), result)


def test_map_aa_code_with_standard_stop_codon():
    # Stop codons are returned unchanged
    input = '*'
    result = '*'
    assert_equals(validate_annotated_vcfs.map_aa_codes(input), result)


def test_map_aa_code_with_non_standard_stop_codon():
    # Alternative stop codon notations are all mapped to *
    input = 'X'
    result = '*'
    assert_equals(validate_annotated_vcfs.map_aa_codes(input), result)


def test_map_aa_code_with_del_codon_notation():
    # del should be returned as is, representing deletion
    input = 'del'
    result = 'DEL'
    assert_equals(validate_annotated_vcfs.map_aa_codes(input), result)


#######################################################################################################################
# Tests for get_unique_tagged_entry_values
#######################################################################################################################
def test_get_tagged_entry_value_with_single_item():
    input = 'LAA_CHANGE=ENTRY'
    tag = 'LAA_CHANGE'
    result = ['ENTRY']
    assert_items_equal(validate_annotated_vcfs.get_unique_tagged_entry_values(input, tag), result)


def test_get_tagged_entry_value_with_multiple_entries():
    input = 'LAA_CHANGE=ENTRY;AA_CHANGE=ENTRY2'
    tag = 'LAA_CHANGE'
    result = ['ENTRY']
    assert_items_equal(validate_annotated_vcfs.get_unique_tagged_entry_values(input, tag), result)


def test_get_tagged_entry_value_with_multiple_entries_get_second_item():
    input = 'LAA_CHANGE=ENTRY;SECOND_TAG=ENTRY2'
    tag = 'SECOND_TAG'
    result = ['ENTRY2']
    assert_items_equal(validate_annotated_vcfs.get_unique_tagged_entry_values(input, tag), result)


def test_get_tagged_entry_value_with_substring_tag():
    # Multiple entries, one is substring of another. Require exact tag match.
    input = 'LAA_CHANGE=ENTRY;AA_ENTRY=ENTRY2'
    tag = 'AA_ENTRY'
    result = ['ENTRY2']
    assert_items_equal(validate_annotated_vcfs.get_unique_tagged_entry_values(input, tag), result)


def test_get_tagged_entry_value_with_nonexisting_tag():
    input = 'LAA_CHANGE=ENTRY;AA_CHANGE=ENTRY2'
    tag = 'NOT_IN_LIST'
    result = []
    assert_raises(ValueError, validate_annotated_vcfs.get_unique_tagged_entry_values, input, tag)


def test_get_tagged_entry_value_with_empty_tag():
    # Request a tag with an empty entry
    input = 'LAA_CHANGE=;AA_ENTRY=ENTRY2'
    tag = 'LAA_CHANGE'
    result = []
    assert_items_equal(validate_annotated_vcfs.get_unique_tagged_entry_values(input, tag), result)


def test_get_tagged_entry_values_with_values_list():
    # Test cases where there is a list of entries as values for in tag-value pairs
    input = 'LAA_CHANGE=A,B,C,D;AA_ENTRY=E,F,G,H'
    tag = 'LAA_CHANGE'
    result = ['A', 'B', 'C', 'D']
    assert_items_equal(validate_annotated_vcfs.get_unique_tagged_entry_values(input, tag), result)


#######################################################################################################################
# Tests for remove_p_dot_notation
#######################################################################################################################
def test_remove_p_dot_notation_with_standard_notation():
    # Basic notation with no parentheses
    input = 'p.Gly47Arg'
    result = 'Gly47Arg'
    assert_equals(validate_annotated_vcfs.remove_p_dot_notation(input), result)


def test_remove_p_dot_notation_with_parentheses():
    # Notation with enclosing parentheses
    input = 'p.(Lys5799Glu)'
    result = 'Lys5799Glu'
    assert_equals(validate_annotated_vcfs.remove_p_dot_notation(input), result)


def test_remove_p_dot_notation_with_brackets():
    # Use of brackets instead of parentheses
    input = 'p.[Lys5799Glu]'
    result = 'Lys5799Glu'
    assert_equals(validate_annotated_vcfs.remove_p_dot_notation(input), result)


def test_remove_p_dot_notation_with_missing_opening_parentheses():
    # Missing starting parentheses
    input = 'p.(Met563Lys'
    result = 'Met563Lys'
    assert_equals(validate_annotated_vcfs.remove_p_dot_notation(input), result)


def test_remove_p_dot_notation_with_missing_closing_parentheses():
    # Missing closing parentheses
    input = 'p.Met563Lys)'
    result = 'Met563Lys'
    assert_equals(validate_annotated_vcfs.remove_p_dot_notation(input), result)


#######################################################################################################################
# Tests for get_laa_change
#######################################################################################################################
def test_get_laa_change_with_standard_input():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'LAA_CHANGE=p.(TYR457GLY)'
    result = ('TYR', 'GLY')
    assert_equals(validate_annotated_vcfs.get_laa_change(input), result)


def test_get_laa_change_with_stop_codon():
    # Test stop codon notation
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'LAA_CHANGE=p.(TYR457*)'
    result = ('TYR', '*')
    assert_equals(validate_annotated_vcfs.get_laa_change(input), result)


def test_get_laa_change_with_alternate_stop_codon():
    # Test alternate stop codon notation (X rather than *)
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'LAA_CHANGE=p.(TYR457X)'
    result = ('TYR', 'X')
    assert_equals(validate_annotated_vcfs.get_laa_change(input), result)


def test_get_laa_change_with_invalid_notation():
    # Should fail with amino acid sequences that are more than 3 letters long (must be invalid)
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'LAA_CHANGE=p.(fs*)'
    assert_raises(ValueError, validate_annotated_vcfs.get_laa_change, input)


#######################################################################################################################
# Tests for get_vep_aa_change
#######################################################################################################################
def test_get_vep_aa_change_with_standard_list():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'AA_CHANGE=W/T,A/T'

    result = [('W', 'T'), ('A', 'T')]
    assert_items_equal(validate_annotated_vcfs.get_vep_aa_change(input), result)


def test_get_vep_aa_change_with_synonymous_mutation():
    input = [''] * (validate_annotated_vcfs.INFO_COLUMN_INDEX + 1)
    input[validate_annotated_vcfs.INFO_COLUMN_INDEX] = 'AA_CHANGE=W/T,A'

    result = [('W', 'T'), ('A', 'A')]
    assert_items_equal(validate_annotated_vcfs.get_vep_aa_change(input), result)


#######################################################################################################################
# Tests for get_26K_allele_frequency
#######################################################################################################################

# TODO - not implemented




