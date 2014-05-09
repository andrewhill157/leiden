from nose.tools import assert_equals
from ..leiden.leiden.leiden_database import *


class TestLeidenDatabase():

    @classmethod
    def setup_class(cls):
        cls.database = LeidenDatabase('http://www.dmd.nl/nmdb2/')

    def test_extract_LOVD_version_number(cls):
        # LOVD2 Installation
        input = 'http://www.dmd.nl/nmdb2/'
        result = 2
        assert_equals(LeidenDatabase.extract_LOVD_version_number(input), result)

        # LOVD3 Installation
        input = 'http://mseqdr.lumc.edu/GEDI/'
        result = 3
        assert_equals(LeidenDatabase.extract_LOVD_version_number(input), result)

    def test_get_version_number(cls):
        # Parent class, does not have a version number
        assert_equals(cls.database.get_version_number(), '')

    def test_set_gene_id(cls):
        # TODO not sure what test would be appropriate here
        pass

    def test_get_variant_database_url(cls):
        # Abstract, not implemented for parent class
        pass

    def test_gene_homepage_url(cls):
        # Abstract, not implemented for parent class
        pass

    def test_get_available_genes(cls):
        # Abstract, not implemented for parent class
        pass

    def test_get_transcript_refseqid(cls):
        # Will be tested separately for each subclass
        pass

    def test_get_table_headers(cls):
        # Abstract, not implemented for parent class
        pass

    def test_get_table_data(cls):
        #TODO not sure what test would be appropriate here
        pass


class TestLOVD2Database():

    @classmethod
    def setup_class(cls):
        # Two LOVD2 installations
        cls.database1 = LOVD2Database('http://www.dmd.nl/nmdb2/')
        cls.database2 = LOVD2Database('http://grenada.lumc.nl/LOVD2/eye/')

    def test_get_lovd_version(cls):
        # Static method, not overridden in subclasses
        pass

    def test_get_version_number(cls):
        # LOVD2 is version 2
        result = 2
        assert_equals(cls.database1.get_version_number(), result)

    def test_set_gene_id(cls):
        # Not overridden in this class
        pass

    def test_get_variant_database_url(cls):
        # Two genes from database1
        input = 'ACTA1'
        result = 'http://www.dmd.nl/nmdb2/variants.php?action=search_unique&select_db=ACTA1&limit=1000'
        cls.database1.set_gene_id(input)
        assert_equals(cls.database1.get_variant_database_url(), result)

        input = 'CAPN3'
        result = 'http://www.dmd.nl/nmdb2/variants.php?action=search_unique&select_db=CAPN3&limit=1000'
        cls.database1.set_gene_id(input)
        assert_equals(cls.database1.get_variant_database_url(), result)

        # Two genes from database2
        input = 'BBS2'
        result = 'http://grenada.lumc.nl/LOVD2/eye/variants.php?action=search_unique&select_db=BBS2&limit=1000'
        cls.database2.set_gene_id(input)
        assert_equals(cls.database2.get_variant_database_url(), result)

        input = 'CRYAA'
        result = 'http://grenada.lumc.nl/LOVD2/eye/variants.php?action=search_unique&select_db=CRYAA&limit=1000'
        cls.database2.set_gene_id(input)
        assert_equals(cls.database2.get_variant_database_url(), result)

    def test_get_gene_homepage_url(cls):
        # Two genes from database1
        input = 'ACTA1'
        result = 'http://www.dmd.nl/nmdb2/home.php?select_db=ACTA1'
        cls.database1.set_gene_id(input)
        assert_equals(cls.database1.get_gene_homepage_url(), result)

        input = 'CAPN3'
        result = 'http://www.dmd.nl/nmdb2/home.php?select_db=CAPN3'
        cls.database1.set_gene_id(input)
        assert_equals(cls.database1.get_gene_homepage_url(), result)

        # Two genes from database2
        input = 'BBS2'
        result = 'http://grenada.lumc.nl/LOVD2/eye/home.php?select_db=BBS2'
        cls.database2.set_gene_id(input)
        assert_equals(cls.database2.get_gene_homepage_url(), result)

        input = 'CRYAA'
        result = 'http://grenada.lumc.nl/LOVD2/eye/home.php?select_db=CRYAA'
        cls.database2.set_gene_id(input)
        assert_equals(cls.database2.get_gene_homepage_url(), result)

    def test_get_available_genes(cls):
        # TODO not sure what test would be appropriate here
        pass

    def test_get_transcript_refseqid(cls):
        # Two genes from database1
        input = 'ACTA1'
        result = 'NM_001100.3'
        cls.database1.set_gene_id(input)
        assert_equals(cls.database1.get_transcript_refseqid(), result)

        input = 'CAPN3'
        result = 'NM_000070.2'
        cls.database1.set_gene_id(input)
        assert_equals(cls.database1.get_transcript_refseqid(), result)

        # Two genes from database2
        input = 'BBS2'
        result = 'NM_031885.3'
        cls.database2.set_gene_id(input)
        assert_equals(cls.database2.get_transcript_refseqid(), result)

        input = 'CRYAA'
        result = 'NM_000394.2'
        cls.database2.set_gene_id(input)
        assert_equals(cls.database2.get_transcript_refseqid(), result)

    def test_get_table_headers(cls):
        #TODO not sure what test would be appropriate here
        pass

    def test_get_table_data_page_n(cls):
        #TODO not sure what test would be appropriate here
        pass


class TestLOVD3Database():

    @classmethod
    def setup_class(cls):
        cls.database1 = LOVD3Database('http://mseqdr.lumc.edu/GEDI/')

    def test_get_lovd_version(cls):
        # Not overridden in this class
        pass

    def test_get_version_number(cls):
        # LOVD3 is version 3
        result = 3
        assert_equals(cls.database1.get_version_number(), result)

    def test_set_gene_id(cls):
        # Not overridden in this class
        pass

    def test_get_variant_database_url(cls):
        # Two genes from database1
        input = 'BBS1'
        result = 'http://mseqdr.lumc.edu/GEDI/variants/BBS1?page_size=1000&page=1'
        cls.database1.set_gene_id(input)
        assert_equals(cls.database1.get_variant_database_url(), result)

        input = 'CTC1'
        result = 'http://mseqdr.lumc.edu/GEDI/variants/CTC1?page_size=1000&page=1'
        cls.database1.set_gene_id(input)
        assert_equals(cls.database1.get_variant_database_url(), result)

    def test_get_gene_homepage_url(cls):
        # Two genes from database1
        input = 'BBS1'
        result = 'http://mseqdr.lumc.edu/GEDI/genes/BBS1?page_size=1000&page=1'
        cls.database1.set_gene_id(input)
        assert_equals(cls.database1.get_gene_homepage_url(), result)

        input = 'CTC1'
        result = 'http://mseqdr.lumc.edu/GEDI/genes/CTC1?page_size=1000&page=1'
        cls.database1.set_gene_id(input)
        assert_equals(cls.database1.get_gene_homepage_url(), result)

    def test_get_available_genes(cls):
        # TODO not sure what test would be appropriate here
        pass

    def test_get_transcript_refseqid(cls):
        # Two genes from database1
        input = 'BBS1'
        result = 'NM_024649.4'
        cls.database1.set_gene_id(input)
        assert_equals(cls.database1.get_transcript_refseqid(), result)

        input = 'CTC1'
        result = 'NM_025099.5'
        cls.database1.set_gene_id(input)
        assert_equals(cls.database1.get_transcript_refseqid(), result)

    def test_get_table_headers(cls):
        #TODO not sure what test would be appropriate here
        pass

    def test_get_table_data_page_n(cls):
        #TODO not sure what test would be appropriate here
        pass
