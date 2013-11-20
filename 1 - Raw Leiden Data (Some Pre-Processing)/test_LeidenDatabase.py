import unittest
from unittest import TestCase
import Leiden


__author__ = 'andrewhill'


class TestLeidenDatabase(TestCase):
    def setUp(self):
        self.database = Leiden.LeidenDatabase("ACTA1")
        self.database2 = Leiden.LeidenDatabase("DYSF")

    def test_get_transcript_refseqid(self):
        refseqid = self.database.get_transcript_refseqid()
        refseqid2 = self.database2.get_transcript_refseqid()
        self.assertEquals(refseqid, "NM_001100.3")
        self.assertEquals(refseqid2, "NM_003494.3")

if __name__ == '__main__':
    unittest.main()