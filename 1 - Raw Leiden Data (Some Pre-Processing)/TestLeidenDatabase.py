import unittest  # Testing framework
from LeidenDatabase import *  # Module under test


class TestTextProcessing(unittest.TestCase):

    def test_get_pmid(self):
        # Typical PubMed URL
        input = 'http://www.ncbi.nlm.nih.gov/pubmed/19562689'
        self.assertEqual(TextProcessing.get_pmid(input), '19562689')

        # Early PubMed URL with shorter ID
        input = 'http://www.ncbi.nlm.nih.gov/pubmed/1592'
        self.assertEqual(TextProcessing.get_pmid(input), '1592')

        # PMID is below the 4-digit minimum, should raise exception
        input = 'http://www.ncbi.nlm.nih.gov/pubmed/34'
        self.assertRaises(ValueError, TextProcessing.get_pmid, input)

        # No PMID present, should raise exception
        input = 'http://www.ncbi.nlm.nih.gov/pubmed/'
        self.assertRaises(ValueError, TextProcessing.get_pmid, input)

    def test_get_omimid(self):
        # Typical OMIM URL
        input = 'http://www.omim.org/entry/102610#0001'
        self.assertEqual(TextProcessing.get_omimid(input), '102610#0001')

        # Shortened, but valid, OMIM URL
        input = 'http://www.omim.org/entry/0#0'
        self.assertEqual(TextProcessing.get_omimid(input), '0#0')

        # Invalid OMIM URL - no digits after #
        input = 'http://www.omim.org/entry/102610#'
        self.assertRaises(ValueError, TextProcessing.get_omimid, input)

        # Invalid OMIM URL - no digits before #
        input = 'http://www.omim.org/entry/#102610'
        self.assertRaises(ValueError, TextProcessing.get_omimid, input)

        # Invalid OMIM URL - no #
        input = 'http://www.omim.org/entry/102610'
        self.assertRaises(ValueError, TextProcessing.get_omimid, input)

        # Invalid OMIM URL - no digits
        input = 'http://www.omim.org/entry/'
        self.assertRaises(ValueError, TextProcessing.get_omimid, input)

    def test_remove_times_reported(self):
        # Basic test case
        input = 'c.5235A>G (Reported 3 Times)'
        self.assertEqual(TextProcessing.remove_times_reported(input), 'c.5235A>G')

        # Basic test case with more digits in times reported
        input = 'c.5235A>G (Reported 403 Times)'
        self.assertEqual(TextProcessing.remove_times_reported(input), 'c.5235A>G')

        # Alternate position for target text
        input = '(Reported 403 Times) c.5235A>G'
        self.assertEqual(TextProcessing.remove_times_reported(input), 'c.5235A>G')

        # Comparison is not case sensitive
        input = 'c.5235A>G (rePoRted 403 times)'
        self.assertEqual(TextProcessing.remove_times_reported(input), 'c.5235A>G')

        # Should return unchanged without altering white-space
        input = ' c.5235A>G '
        self.assertEqual(TextProcessing.remove_times_reported(input), ' c.5235A>G ')

    def test_find_string_index(self):
        # Basic test, should return first instance of strings appearing twice
        input = ['test', 'other', 'next', 'unit', 'other']
        self.assertEqual(TextProcessing.find_string_index(input, 'other'), 1)

        # Comparisons should not be case or white-space sensitive
        input = ['Test ', 'OtHeR ', ' nExt ', 'unit']
        self.assertEqual(TextProcessing.find_string_index(input, 'other'), 1)
        self.assertEqual(TextProcessing.find_string_index(input, 'next '), 2)

        # Comparisons should find substrings
        input = ['Word now', 'Word next', 'next', 'unit']
        self.assertEqual(TextProcessing.find_string_index(input, 'Word'), 0)

        # Return -1 if search string is not found
        input = ['test', 'other', 'next', 'unit', 'other']
        self.assertEqual(TextProcessing.find_string_index(input, 'not in list'), -1)

        # Return -1 if list is empty
        input = []
        self.assertEqual(TextProcessing.find_string_index(input, 'target'), -1)


class TestVariantRemapper(unittest.TestCase):
    def setUp(self):
        self.remapper = VariantRemapper()

    def test_remap_variant(self):
        # Simple SNP remapping example
        input = 'NM_001100.3:c.24C>A'
        result = 'NC_000001.10:g.229568839G>T'
        self.assertEqual(self.remapper.remap_variant(input), result)

        # More complex remapping example
        input = 'NM_001100.3:c.-66_-65delinsTC'
        result = 'NC_000001.10:g.229569803_229569804delinsGA'
        self.assertEqual(self.remapper.remap_variant(input), result)

        # Should result in error, empty string returned
        input = 'NM_001100.3:c.='
        result = ''
        self.assertEqual(self.remapper.remap_variant(input), result)

    def test_submit_variant_batch(self):
        self.fail('Not implemented yet.')

    def test_entries_remaining_in_batch(self):
        self.fail('Not implemented yet.')

    def test_get_batch_results(self):
        self.fail('Not implemented yet.')


class TestLeidenDatabase(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_lovd_version(self):
        self.fail('Not implemented yet.')

    def test_get_version_number(self):
        self.fail('Not implemented yet.')

    def test_set_gene_id(self):
        self.fail('Not implemented yet.')

    def test_get_variant_database_url(self):
        self.fail('Not implemented yet.')

    def test_gene_homepage_url(self):
        self.fail('Not implemented yet.')

    def test_get_available_genes(self):
        self.fail('Not implemented yet.')

    def test_get_link_info(self):
        self.fail('Not implemented yet.')

    def test_get_transcript_refseqid(self):
        self.fail('Not implemented yet.')

    def test_get_table_headers(self):
        self.fail('Not implemented yet.')

    def test_get_table_data(self):
        self.fail('Not implemented yet.')

    def test_get_gene_name(self):
        self.fail('Not implemented yet.')


class TestLOVD2Database(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_lovd_version(self):
        self.fail('Not implemented yet.')

    def test_get_version_number(self):
        self.fail('Not implemented yet.')

    def test_set_gene_id(self):
        self.fail('Not implemented yet.')

    def test_get_variant_database_url(self):
        self.fail('Not implemented yet.')

    def test_gene_homepage_url(self):
        self.fail('Not implemented yet.')

    def test_get_available_genes(self):
        self.fail('Not implemented yet.')

    def test_get_link_info(self):
        self.fail('Not implemented yet.')

    def test_get_transcript_refseqid(self):
        self.fail('Not implemented yet.')

    def test_get_table_headers(self):
        self.fail('Not implemented yet.')

    def test_get_table_data(self):
        self.fail('Not implemented yet.')

    def test_get_gene_name(self):
        self.fail('Not implemented yet.')


class TestLOVD3Database(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_lovd_version(self):
        self.fail('Not implemented yet.')

    def test_get_version_number(self):
        self.fail('Not implemented yet.')

    def test_set_gene_id(self):
        self.fail('Not implemented yet.')

    def test_get_variant_database_url(self):
        self.fail('Not implemented yet.')

    def test_gene_homepage_url(self):
        self.fail('Not implemented yet.')

    def test_get_available_genes(self):
        self.fail('Not implemented yet.')

    def test_get_link_info(self):
        self.fail('Not implemented yet.')

    def test_get_transcript_refseqid(self):
        self.fail('Not implemented yet.')

    def test_get_table_headers(self):
        self.fail('Not implemented yet.')

    def test_get_table_data(self):
        self.fail('Not implemented yet.')

    def test_get_gene_name(self):
        self.fail('Not implemented yet.')

if __name__ == '__main__':
    unittest.main()