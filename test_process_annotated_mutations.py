from nose.tools import assert_equals
from nose.tools import assert_raises
from nose.tools import assert_items_equal
import process_annotated_mutations


#######################################################################################################################
# Tests for has_vep_aa_change
#######################################################################################################################


#######################################################################################################################
# Tests for has_lovd_aa_change
#######################################################################################################################


#######################################################################################################################
# Tests for get_severe_impact
#######################################################################################################################


#######################################################################################################################
# Tests for get_severe_impact
#######################################################################################################################


#######################################################################################################################
# Tests for is_concordant
#######################################################################################################################


#######################################################################################################################
# Tests for is_concordant_splice_mutation
#######################################################################################################################


#######################################################################################################################
# Tests for is_concordant_inframe_codon_loss
#######################################################################################################################


#######################################################################################################################
# Tests for plot_allele_frequency_histogram
#######################################################################################################################


#######################################################################################################################
# Tests for get_ucsc_location_link
#######################################################################################################################
def test_get_ucsc_location_link_with_standard_interval():
    # Basic test to ensure URL is constructed correctly
    chromosome = '1'
    start_coordinate = '15613'
    end_coordinate = '52552'
    result = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1%3A15613-52552'

    assert_equals(process_annotated_mutations.get_ucsc_location_link(chromosome, start_coordinate, end_coordinate), result)


#######################################################################################################################
# Tests for map_aa_codes
#######################################################################################################################
def test_map_aa_codes_with_three_letter_code():
    # Remap three letter code to one letter code
    input = 'GLU'
    result = 'E'
    assert_equals(process_annotated_mutations.map_aa_codes(input), result)


def test_mapp_aa_code_with_lower_case():
    # Case should not matter
    input = 'met'
    result = 'M'
    assert_equals(process_annotated_mutations.map_aa_codes(input), result)


def test_mapp_aa_code_with_one_letter_code():
    # One letter codes returned unchanged
    input = 'E'
    result = 'E'
    assert_equals(process_annotated_mutations.map_aa_codes(input), result)


def test_mapp_aa_code_with_lower_case_one_letter_code():
    # One letter codes are returned as capitals regardless of input case
    input = 'k'
    result = 'K'
    assert_equals(process_annotated_mutations.map_aa_codes(input), result)


def test_mapp_aa_code_with_standard_stop_codon():
    # Stop codons are returned unchanged
    input = '*'
    result = '*'
    assert_equals(process_annotated_mutations.map_aa_codes(input), result)


def test_map_aa_code_with_non_standard_stop_codon():
    # Alternative stop codon notations are all mapped to *
    input = 'X'
    result = '*'
    assert_equals(process_annotated_mutations.map_aa_codes(input), result)


def test_mapp_aa_code_with_del_codon_notation():
    # del should be returned as is, representing deletion
    input = 'del'
    result = 'DEL'
    assert_equals(process_annotated_mutations.map_aa_codes(input), result)


#######################################################################################################################
# Tests for get_unique_tagged_entry_values
#######################################################################################################################
def test_get_tagged_entry_value_with_single_item():
    input = 'LAA_CHANGE=ENTRY'
    tag = 'LAA_CHANGE'
    result = ['ENTRY']
    assert_items_equal(process_annotated_mutations.get_unique_tagged_entry_values(input, tag), result)


def test_get_tagged_entry_value_with_multiple_entries():
    input = 'LAA_CHANGE=ENTRY;AA_CHANGE=ENTRY2'
    tag = 'LAA_CHANGE'
    result = ['ENTRY']
    assert_items_equal(process_annotated_mutations.get_unique_tagged_entry_values(input, tag), result)


def test_get_tagged_entry_value_with_multiple_entries_get_second_item():
    input = 'LAA_CHANGE=ENTRY;SECOND_TAG=ENTRY2'
    tag = 'SECOND_TAG'
    result = ['ENTRY2']
    assert_items_equal(process_annotated_mutations.get_unique_tagged_entry_values(input, tag), result)


def test_get_tagged_entry_value_with_substring_tag():
    # Multiple entries, one is substring of another. Require exact tag match.
    input = 'LAA_CHANGE=ENTRY;AA_ENTRY=ENTRY2'
    tag = 'AA_ENTRY'
    result = ['ENTRY2']
    assert_items_equal(process_annotated_mutations.get_unique_tagged_entry_values(input, tag), result)


def test_get_tagged_entry_value_with_nonexisting_tag():
    input = 'LAA_CHANGE=ENTRY;AA_CHANGE=ENTRY2'
    tag = 'NOT_IN_LIST'
    result = []
    assert_raises(ValueError, process_annotated_mutations.get_unique_tagged_entry_values, input, tag)


def test_get_tagged_entry_value_with_empty_tag():
    # Request a tag with an empty entry
    input = 'LAA_CHANGE=;AA_ENTRY=ENTRY2'
    tag = 'LAA_CHANGE'
    result = []
    assert_items_equal(process_annotated_mutations.get_unique_tagged_entry_values(input, tag), result)


def test_get_tagged_entry_values_with_values_list():
    # Test cases where there is a list of entries as values for in tag-value pairs
    input = 'LAA_CHANGE=A,B,C,D;AA_ENTRY=E,F,G,H'
    tag = 'LAA_CHANGE'
    result = ['A', 'B', 'C', 'D']
    assert_items_equal(process_annotated_mutations.get_unique_tagged_entry_values(input, tag), result)


#######################################################################################################################
# Tests for remove_p_dot_notation
#######################################################################################################################
def test_remove_p_dot_notation_with_standard_notation():
    # Basic notation with no parentheses
    input = 'p.Gly47Arg'
    result = 'Gly47Arg'
    assert_equals(process_annotated_mutations.remove_p_dot_notation(input), result)


def test_remove_p_dot_notation_with_parentheses():
    # Notation with enclosing parentheses
    input = 'p.(Lys5799Glu)'
    result = 'Lys5799Glu'
    assert_equals(process_annotated_mutations.remove_p_dot_notation(input), result)


def test_remove_p_dot_notation_with_brackets():
    # Use of brackets instead of parentheses
    input = 'p.[Lys5799Glu]'
    result = 'Lys5799Glu'
    assert_equals(process_annotated_mutations.remove_p_dot_notation(input), result)


def test_remove_p_dot_notation_with_missing_opening_parentheses():
    # Missing starting parentheses
    input = 'p.(Met563Lys'
    result = 'Met563Lys'
    assert_equals(process_annotated_mutations.remove_p_dot_notation(input), result)


def test_remove_p_dot_notation_with_missing_closing_parentheses():
    # Missing closing parentheses
    input = 'p.Met563Lys)'
    result = 'Met563Lys'
    assert_equals(process_annotated_mutations.remove_p_dot_notation(input), result)


#######################################################################################################################
# Tests for get_laa_change
#######################################################################################################################


#######################################################################################################################
# Tests for get_vep_aa_change
#######################################################################################################################
def test_get_vep_aa_change_with_standard_list():
    input = [''] * (process_annotated_mutations.INFO_COLUMN_INDEX + 1)
    input[process_annotated_mutations.INFO_COLUMN_INDEX] = 'AA_CHANGE=W/T,A/T'

    result = [('W', 'T'), ('A', 'T')]
    assert_items_equal(process_annotated_mutations.get_vep_aa_change(input), result)


def test_get_vep_aa_change_with_synonymous_mutation():
    input = [''] * (process_annotated_mutations.INFO_COLUMN_INDEX + 1)
    input[process_annotated_mutations.INFO_COLUMN_INDEX] = 'AA_CHANGE=W/T,A'

    result = [('W', 'T'), ('A', 'A')]
    assert_items_equal(process_annotated_mutations.get_vep_aa_change(input), result)


#######################################################################################################################
# Tests for get_26K_allele_frequency
#######################################################################################################################





