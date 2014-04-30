"""
Andrew Hill
MacArthur Lab - 2014

Functions for handling remapping between HGVS and genomic coordinates.
"""

import os
import re

import hgvs
import hgvs.utils
from pygr.seqdb import SequenceFileDB

from leiden.input_output import file_io
from leiden.lovd import utilities


class VariantRemapper:
    """
    Class containing functions for remapping of variants from HGVS to genomic coordinate notation.
    """

    def __init__(self):
        """
        Initializes hg19 reference and reference transcripts
        """

        genome_path = os.path.join(os.path.dirname(__file__), 'resources', 'hg19.fa')
        refseq_path = os.path.join(os.path.dirname(__file__), 'resources', 'genes.refGene')

        # Read genome sequence using pygr.
        self.genome = SequenceFileDB(genome_path)

        # Read RefSeq transcripts into a python dict.
        with open(refseq_path) as infile:
            self.transcripts = hgvs.utils.read_transcripts(infile)

    def hgvs_to_vcf(self, hgvs_variant):
        """
        Converts a single variant provided in HGVS notation to genomic coordinate notation.

        See U(https://humgenprojects.lumc.nl/trac/mutalyzer/wiki/PositionConverter) for more information on acceptable
        inputs and outputs, and remapping functionality.

        @param hgvs_variant: HGVS description of variant, such as NM_001100.3:c.137T>C. The portion prior to the colon is
        the refseqID used as the reference for the variant The portion after the colon is an HGVS-style description
        of the mutation (a SNP from T to C at location 137 in the example above.
        @type variant: string
        @return: A tuple (chromosome_number, coordinate, ref, alt) in that order denoting the VCF notation of the variant
        @rtype: tuple of strings
        """

        # Library requires string not unicode, ensure format is correct
        hgvs_variant = str(hgvs_variant)

        chromosome_number, coordinate, ref, alt = hgvs.parse_hgvs_name(hgvs_variant, self.genome, get_transcript=self._get_transcript)
        chromosome_number = re.match('chr(.+)', chromosome_number).group(1)
        coordinate = str(coordinate)

        return chromosome_number, coordinate, ref, alt

    def vcf_to_hgvs(self, reference_transcript, vcf_notation):
        """
        Converts a single VCF notation variant to HGVS notation relative to a given transcript.

        @param reference_transcript: the refseq id of the reference transcript to use for HGVS notation
        @type reference_transcript: string
        @param vcf_notation: a tuple containing elements chromosome_number, coordinate, ref, and alt in that order
        @type vcf_notation: tuple of strings
        @return: hgvs notatation of variant in format reference_transcript:hgvs_description
        @rtype: string
        """

        chromosome_number, coordinate, ref, alt = vcf_notation
        coordinate = int(coordinate)

        transcript = self._get_transcript(reference_transcript)

        return hgvs.format_hgvs_name(chromosome_number, coordinate, ref, alt, self.genome, transcript)

    def _get_transcript(self, name):
        """
        Callback to provide reference transcript by its name

        @param name: name of reference transcript
        @type name: string
        @return: line of information on transcript from resource file
        """
        return self.transcripts.get(name)


def generate_vcf_from_hgvs(input_file, output_file, hgvs_column, protein_change_column, column_delimiter='\t'):
    """
    Generate VCF files from files containing variants in HGVS notation. First row must contain column labels.

    @param input_file: path to input file containing HGVS variants.
    @type input_file: string
    @param output_file: path to output file for VCF file
    @type output_file: string
    @param hgvs_column: column label for HGVS notation column
    @type hgvs_column: string
    @param protein_change_column: column label for predicted protein change column
    @type protein_change_column: string
    @param column_delimiter: column delimiter to use for input file
    @type column_delimiter: string
    @return: list of remapping error information. Each entry is [file, hgvs_notation, error].
    @rtype: 2D list
    """

    remapper = VariantRemapper()
    remapping_errors = []

    table_data = file_io.read_table_from_file(input_file, column_delimiter=column_delimiter)
    header = table_data.pop(0)

    # Isolate data columns with HGVS mutations and protein change
    hgvs_notation_index = utilities.find_string_index(header, hgvs_column)
    hgvs_notation = []

    protein_change_index = utilities.find_string_index(header, protein_change_column)
    protein_change = []

    # Remap Variants and build list for INFO column tags
    vcf_notation_variants = []

    for row in table_data:
        try:
            vcf_notation = remapper.hgvs_to_vcf(row[hgvs_notation_index])
            vcf_notation_variants.append(vcf_notation)

            hgvs_notation.append(row[hgvs_notation_index])
            protein_change.append(row[protein_change_index])

        except Exception as e:
            remapping_errors.append([file, row[hgvs_notation_index], str(e)])

    info_column_tags = {'HGVS': ('string', 'LOVD HGVS notation describing DNA change', hgvs_notation),
            'LAA_CHANGE': ('string', 'LOVD amino acid change', protein_change)}

    # Write VCF file
    vcf_file_name = os.path.splitext(file)[0] + '.vcf'

    vcf_file_text = file_io.format_vcf_text(vcf_notation_variants, info_column_tags)
    file_io.write_table_to_file(vcf_file_name, vcf_file_text)

    return remapping_errors