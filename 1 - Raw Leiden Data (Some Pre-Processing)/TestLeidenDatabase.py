import unittest  # Testing framework
from LeidenDatabase import *  # Module under test


class TestTextProcessing(unittest.TestCase):

    def test_get_pmid(self):
        # Typical PubMed URL
        input = 'http://www.ncbi.nlm.nih.gov/pubmed/19562689'
        result = '19562689'
        self.assertEqual(TextProcessing.get_pmid(input), result)

        # Early PubMed URL with shorter ID
        input = 'http://www.ncbi.nlm.nih.gov/pubmed/1592'
        result = '1592'
        self.assertEqual(TextProcessing.get_pmid(input), result)

        # PMID is below the 4-digit minimum, should raise exception
        input = 'http://www.ncbi.nlm.nih.gov/pubmed/34'
        self.assertRaises(ValueError, TextProcessing.get_pmid, input)

        # No PMID present, should raise exception
        input = 'http://www.ncbi.nlm.nih.gov/pubmed/'
        self.assertRaises(ValueError, TextProcessing.get_pmid, input)

    def test_get_omimid(self):
        # Typical OMIM URL
        input = 'http://www.omim.org/entry/102610#0001'
        result = '102610#0001'
        self.assertEqual(TextProcessing.get_omimid(input), result)

        # Shortened, but valid, OMIM URL
        input = 'http://www.omim.org/entry/0#0'
        result = '0#0'
        self.assertEqual(TextProcessing.get_omimid(input), result)

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
        result = 'c.5235A>G'
        self.assertEqual(TextProcessing.remove_times_reported(input), result)

        # Basic test case with more digits in times reported
        input = 'c.5235A>G (Reported 403 Times)'
        result = 'c.5235A>G'
        self.assertEqual(TextProcessing.remove_times_reported(input), result)

        # Alternate position for target text
        input = '(Reported 403 Times) c.5235A>G'
        result = 'c.5235A>G'
        self.assertEqual(TextProcessing.remove_times_reported(input), result)

        # Comparison is not case sensitive
        input = 'c.5235A>G (rePoRted 403 times)'
        result = 'c.5235A>G'
        self.assertEqual(TextProcessing.remove_times_reported(input), result)

        # Should return unchanged without altering white-space
        input = ' c.5235A>G '
        result = ' c.5235A>G '
        self.assertEqual(TextProcessing.remove_times_reported(input), result)

    def test_find_string_index(self):
        # Basic test, should return first instance of strings appearing twice
        input = ['test', 'other', 'next', 'unit', 'other']
        result = 1
        self.assertEqual(TextProcessing.find_string_index(input, 'other'), result)

        # Comparisons should not be case or white-space sensitive
        input = ['Test ', 'OtHeR ', ' nExt ', 'unit']
        result = 1
        self.assertEqual(TextProcessing.find_string_index(input, 'other'), result)

        result = 2
        self.assertEqual(TextProcessing.find_string_index(input, 'next '), result)

        # Comparisons should find substrings
        input = ['Word now', 'Word next', 'next', 'unit']
        result = 0
        self.assertEqual(TextProcessing.find_string_index(input, 'Word'), result)

        # Return -1 if search string is not found
        input = ['test', 'other', 'next', 'unit', 'other']
        result = -1
        self.assertEqual(TextProcessing.find_string_index(input, 'not in list'), result)

        # Return -1 if list is empty
        input = []
        result = -1
        self.assertEqual(TextProcessing.find_string_index(input, 'target'), result)


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

    def test_get_chromosome_number(self):
        # One digit chromosome number
        input = 'NC_0000001.2:g.524264A>C'
        result = '1'
        self.assertEqual(VariantRemapper.get_chromosome_number(input), result)

        # Two digit chromosome number
        input = 'NC_0000012.1:g.524264A>C'
        result = '12'
        self.assertEqual(VariantRemapper.get_chromosome_number(input), result)

    def test_get_coordinates(self):
        # Long coordinate number
        input = 'NC_0000001.2:g.524264A>C'
        result = '524264'
        self.assertEqual(VariantRemapper.get_coordinates(input), result)

        # Single-digit coordinate number
        input = 'NC_0000012.1:g.1G>T'
        result = '1'
        self.assertEqual(VariantRemapper.get_coordinates(input), result)

    def test_get_ref(self):
        # Example 1
        input = 'NC_0000001.2:g.524264A>C'
        result = 'A'
        self.assertEqual(VariantRemapper.get_ref(input), result)

        # Example 2
        input = 'NC_0000012.1:g.1G>T'
        result = 'G'
        self.assertEqual(VariantRemapper.get_ref(input), result)

    def test_get_alt(self):
        # Example 1
        input = 'NC_0000001.2:g.524264A>C'
        result = 'C'
        self.assertEqual(VariantRemapper.get_alt(input), result)

        # Example 2
        input = 'NC_0000012.1:g.1G>T'
        result = 'T'
        self.assertEqual(VariantRemapper.get_alt(input), result)


class TestLeidenDatabase(unittest.TestCase):
    def setUp(self):
        self.database = LeidenDatabase('http://www.dmd.nl/nmdb2/')

    def test_get_lovd_version(self):
        # LOVD2 Installation
        input = 'http://www.dmd.nl/nmdb2/'
        result = 2
        self.assertEqual(LeidenDatabase.get_lovd_version(input), result)

        # LOVD3 Installation
        input = 'http://mseqdr.lumc.edu/GEDI/'
        result = 3
        self.assertEqual(LeidenDatabase.get_lovd_version(input), result)

    def test_get_version_number(self):
        # Parent class, does not have a version number
        self.assertEqual(self.database.get_version_number(), '')

    def test_set_gene_id(self):
        self.fail('Not implemented yet.')

    def test_get_variant_database_url(self):
        # Abstract, not implemented for parent class
        pass

    def test_gene_homepage_url(self):
        # Abstract, not implemented for parent class
        pass

    def test_get_available_genes(self):
        # Abstract, not implemented for parent class
        pass

    def test_get_link_info(self):
        self.fail('Not implemented yet.')

    def test_get_transcript_refseqid(self):
        # Will be tested separately for each subclass
        pass

    def test_get_table_headers(self):
        # Abstract, not implemented for parent class
        pass

    def test_get_table_data(self):
        # Abstract, not implemented for parent class
        pass


class TestLOVD2Database(unittest.TestCase):
    def setUp(self):
        # Two LOVD2 installations
        self.database1 = LOVD2Database('http://www.dmd.nl/nmdb2/')
        self.database2 = LOVD2Database('http://grenada.lumc.nl/LOVD2/eye/')

    def test_get_lovd_version(self):
        # Static method, not overridden in subclasses
        pass

    def test_get_version_number(self):
        # LOVD2 is version 2
        result = 2
        self.assertEqual(self.database1.get_version_number(), result)

    def test_set_gene_id(self):
        self.fail('Not implemented yet.')

    def test_get_variant_database_url(self):
        # Two genes from database1
        input = 'ACTA1'
        result = 'http://www.dmd.nl/nmdb2/variants.php?action=search_unique&select_db=ACTA1&limit=1000'
        self.assertEqual(self.database1.get_variant_database_url(input), result)

        input = 'CAPN3'
        result = 'http://www.dmd.nl/nmdb2/variants.php?action=search_unique&select_db=CAPN3&limit=1000'
        self.assertEqual(self.database1.get_variant_database_url(input), result)

        # Two genes from database2
        input = 'BBS2'
        result = 'http://grenada.lumc.nl/LOVD2/eye/variants.php?action=search_unique&select_db=BBS2&limit=1000'
        self.assertEqual(self.database2.get_variant_database_url(input), result)

        input = 'CRYAA'
        result = 'http://grenada.lumc.nl/LOVD2/eye/variants.php?action=search_unique&select_db=CRYAA&limit=1000'
        self.assertEqual(self.database2.get_variant_database_url(input), result)

    def test_get_gene_homepage_url(self):
        # Two genes from database1
        input = 'ACTA1'
        result = 'http://www.dmd.nl/nmdb2/home.php?select_db=ACTA1'
        self.assertEqual(self.database1.get_gene_homepage_url(input), result)

        input = 'CAPN3'
        result = 'http://www.dmd.nl/nmdb2/home.php?select_db=CAPN3'
        self.assertEqual(self.database1.get_gene_homepage_url(input), result)

        # Two genes from database2
        input = 'BBS2'
        result = 'http://grenada.lumc.nl/LOVD2/eye/home.php?select_db=BBS2'
        self.assertEqual(self.database2.get_gene_homepage_url(input), result)

        input = 'CRYAA'
        result = 'http://grenada.lumc.nl/LOVD2/eye/home.php?select_db=CRYAA'
        self.assertEqual(self.database2.get_gene_homepage_url(input), result)

    def test_get_available_genes(self):
        self.fail('Not implemented yet.')

    def test_get_link_info(self):
        self.fail('Not implemented yet.')

    def test_get_transcript_refseqid(self):
        # Two genes from database1
        input = 'ACTA1'
        result = 'NM_001100.3'
        self.assertEqual(self.database1.get_transcript_refseqid(input), result)

        input = 'CAPN3'
        result = 'NM_000070.2'
        self.assertEqual(self.database1.get_transcript_refseqid(input), result)

        # Two genes from database2
        input = 'BBS2'
        result = 'NM_031885.3'
        self.assertEqual(self.database2.get_transcript_refseqid(input), result)

        input = 'CRYAA'
        result = 'NM_000394.2'
        self.assertEqual(self.database2.get_transcript_refseqid(input), result)

    def test_get_table_headers(self):
        self.fail('Not implemented yet.')

    def test_get_table_data(self):
        self.fail('Not implemented yet.')


class TestLOVD3Database(unittest.TestCase):
    def setUp(self):
        self.database1 = LOVD3Database('http://mseqdr.lumc.edu/GEDI/')
        self.database2 = LOVD3Database('http://mseqdr.lumc.edu/MITO/')

    def test_get_lovd_version(self):
        # Static method, not overridden in subclasses
        pass

    def test_get_version_number(self):
        # LOVD3 is version 3
        result = 3
        self.assertEqual(self.database1.get_version_number(), result)

    def test_set_gene_id(self):
        self.fail('Not implemented yet.')

    def test_get_variant_database_url(self):
        # Two genes from database1
        input = 'BBS1'
        result = 'http://mseqdr.lumc.edu/GEDI/variants/BBS1?page_size=1000&page=1'
        self.assertEqual(self.database1.get_variant_database_url(input), result)

        input = 'CTC1'
        result = 'http://mseqdr.lumc.edu/GEDI/variants/CTC1?page_size=1000&page=1'
        self.assertEqual(self.database1.get_variant_database_url(input), result)

        # Two genes from database2
        input = 'AMACR'
        result = 'http://mseqdr.lumc.edu/MITO/variants/AMACR?page_size=1000&page=1'
        self.assertEqual(self.database2.get_variant_database_url(input), result)

        input = 'ALDH3A2'
        result = 'http://mseqdr.lumc.edu/MITO/variants/ALDH3A2?page_size=1000&page=1'
        self.assertEqual(self.database2.get_variant_database_url(input), result)

    def test_get_gene_homepage_url(self):
        # Two genes from database1
        input = 'BBS1'
        result = 'http://mseqdr.lumc.edu/GEDI/genes/BBS1?page_size=1000&page=1'
        self.assertEqual(self.database1.get_gene_homepage_url(input), result)

        input = 'CTC1'
        result = 'http://mseqdr.lumc.edu/GEDI/genes/CTC1?page_size=1000&page=1'
        self.assertEqual(self.database1.get_gene_homepage_url(input), result)

        # Two genes from database2
        input = 'AMACR'
        result = 'http://mseqdr.lumc.edu/MITO/genes/AMACR?page_size=1000&page=1'
        self.assertEqual(self.database2.get_gene_homepage_url(input), result)

        input = 'ALDH3A2'
        result = 'http://mseqdr.lumc.edu/MITO/genes/ALDH3A2?page_size=1000&page=1'
        self.assertEqual(self.database2.get_gene_homepage_url(input), result)

    def test_get_available_genes(self):
        self.fail('Not implemented yet.')

    def test_get_link_info(self):
        # Static method, not overridden in subclasses
        pass

    def test_get_transcript_refseqid(self):
        # Two genes from database1
        input = 'BBS1'
        result = 'NM_024649.4'
        self.assertEqual(self.database1.get_transcript_refseqid(input), result)

        input = 'CTC1'
        result = 'NM_025099.5'
        self.assertEqual(self.database1.get_transcript_refseqid(input), result)

        # Two genes from database2
        input = 'AMACR'
        result = 'NM_001167595.1'
        self.assertEqual(self.database2.get_transcript_refseqid(input), result)

        input = 'ALDH3A2'
        result = 'NM_000382.2'
        self.assertEqual(self.database2.get_transcript_refseqid(input), result)

    def test_get_table_headers(self):
        self.fail('Not implemented yet.')

    def test_get_table_data(self):
        self.fail('Not implemented yet.')

if __name__ == '__main__':
    unittest.main()