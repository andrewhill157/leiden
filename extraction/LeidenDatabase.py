from bs4 import BeautifulSoup
import urllib.request
import re
from suds.client import Client
import suds
import base64
import math
from collections import namedtuple


def make_leiden_database(leiden_url):
    """
    Factory method that returns appropriate LeidenDatabase object for LOVD version installed at specified URL.
    Only LOVD2 and LOVD3 installations are supported.

    @param leiden_url: the base URL of the particular Leiden database to be used. For example, the Leiden muscular
    dystrophy pages LOVD2 homepage is http://www.dmd.nl/nmdb2/. This must be a valid URL to base page of database.
    For LOVD3 installations, such as the Genetic Eye Disorder (GEI) Variation Database, the base url will be
    similar to U(http://mseqdr.lumc.edu/GEDI/). Extensions of this URL, such as U(http://mseqdr.lumc.edu/GEDI/genes)
    or U(http://mseqdr.lumc.edu/GEDI/variants) should not be used.
    @type leiden_url: string
    @return: LeidenDatabase object for specified URL.
    @rtype: LeidenDatabase
    @raise: Exception if the LOVD version installed at specified URL is unsupported (anything other than version 2 or 3)
    """

    version = LeidenDatabase.extract_LOVD_version_number(leiden_url)

    # Generate instance of appropriate _LeidenDatabase subclass for installed version
    if version == 2:
        database = LOVD2Database(leiden_url)
    elif version == 3:
        database = LOVD3Database(leiden_url)
    else:
        raise Exception("Unrecognized version number: " + str(version) + "!")

    return database


class Utilities:
    """
    Class containing functions useful for processing of data contained in LOVD installations.\
    All functions are static, class used for organizational purposes only.
    """

    @staticmethod
    def get_pmid(link_url):
        """
        Given a URL to a publication listed on PUBMED, return a string containing the PUBMED ID of the publication.

        @param link_url: URL to the publication on PUBMED. Assumed to be a valid link to a publication on PUBMED.
        For example, U(http://www.ncbi.nlm.nih.gov/pubmed/19562689) is a valid pubmed publication URL. The url must
        contain the PMID in the URL (19562689 in the example here) and contain no other 4 digit or longer numbers.
        @type link_url: string
        @return: PUBMED ID associated with link_url (as specified by the N digit ID included in PUBMED URLs). Assumes
        that PMIDs are at least 4 digits long and that no other 4 digit or longer numeric sequences are contained in
        link_url. Returns empty string if no 4 digit or longer sequence is present in URL.
        @rtype: string
        @raise: ValueError if there is no 4+ digit number in link_url
        """

        # Search for sequences of digits that are four digits or longer in length.
        m = re.compile('\d{4,}')
        results = m.search(link_url)

        # Return entire matched sequence (PMID)
        if results is not None:
            return results.group()
        else:
            raise ValueError('Input URL did not contain 4+ digit number.')

    @staticmethod
    def get_omimid(link_url):
        """
        Given a URL to an entry on OMIM, return a string containing the OMIM ID for the entry.

        @param link_url: URL to the entry on OMIM. Assumed to be a valid link to an entry on PUBMED.
        For example, U(http://www.omim.org/entry/102610#0003) is a valid link to an OMIM entry on the ACTA1 gene.
        The url must contain the gene ID followed by the entry number in the URL separated by a hash mark
        (such as, 102610#0003 in the example URL). URL may not contain other instances of this pattern.
        @type link_url: string
        @return: OMIM entry associated with the URL. This consists of the gene ID (such as 102610 for ACTA1 and a
        specific entry number (0003) separated by a hash mark (102610#0003 in the example above).
        @rtype: string
        @raise: Value Error if link does not contain a valid OMIM ID.
        """

        # Search for sequences of digits separated only by a hash mark
        m = re.compile('\d+#\d+')
        results = m.search(link_url)

        # Return entire matched sequence (OMIM ID)
        if results is not None:
            return results.group()
        else:
            raise ValueError('Input URL did not contain valid OMIM ID.')

    @staticmethod
    def remove_times_reported(hgvs_notation):
        """
        If the hgvs_notation string contains '(Reported N times)' embedded in the HGVS variant description, returns
        a new string with (Reported N times) removed. Comparison is not case sensitive.

        @param hgvs_notation: String, typically an entry in the DNA Change column in table_data for a given variant on
        an LOVD installation.
        @type hgvs_notation: string
        @return: hgvs_notation with instances of (Reported N times) removed. Whitespace surrounding this substring is
        removed in returned string. Return value is equal to input parameter if '(Reported N times)' substring not found.
        @rtype: string
        """
        # Compile case-insensitive regex to match pattern
        m = re.compile('\s*\(Reported \d+ Times\)\s*', re.IGNORECASE)

        # Replace pattern in original string with empty string
        return m.sub('', hgvs_notation)

    @staticmethod
    def find_string_index(string_list, search_string):
        """
        Given a list of strings and a string to search for, returns the index of the first element in the list
        that contains the search string. Note that the comparison is not sensitive to case or leading or trailing
        whitespace characters.

        @param string_list: list of strings
        @type string_list: list of strings
        @param search_string: a string to search for in elements of string_list
        @type search_string: string
        @return: index of the first instance of search_string as a substring of element in string_list. Returns -1 if
        search_string is not found in the string_list.
        @rtype: number
        """

        # Remove leading/trailing whitespace and convert to lowercase before comparisons
        string_list = [x.lower().strip() for x in string_list]
        search_string = search_string.lower().strip()

        for i in range(0, len(string_list)):
            entry = string_list[i]

            if search_string in entry:
                return i

        # search_string not found, return -1
        return -1

    @staticmethod
    def swap(list, i, j):
        """
        Swaps elements at indices i and j in list. Indices must be within bounds of array. Index swapped with itself
        leaves the list unchanged.

        @param list: list of elements
        @type list: list
        @param i: index of element to be swapped with element at index j. Must be within bounds of array.
        @param j: index of element to be swapped with element at index i. Must be within bounds of array.
        @return: list with elements at indices i and j swapped. If i and j are equal, list is unchanged.
        @rtype: list
        """

        list[i], list[j] = list[j], list[i]
        return list

    @staticmethod
    def get_page_html(page_url):
        """
        Returns the html describing the page at the specified URL.

        @param page_url: URL to a specified website
        @type page_url: string
        @return: HTML describing the specified page
        @rtype: string
        @raise: IOError if URL requested website could not be reached
        """
        try:
            with urllib.request.urlopen(page_url) as url:
                url_text = url.read().decode('utf-8', 'ignore')

            return url_text
        except IOError:
            raise IOError('HTML from requested URL could not be retrieved.')


    def deep_copy(nested_list):
        """
        Makes a deep copy of lists that may or may not contain nested lists. Nested items that are not lists will not
        be deep copied, they will be shallow copied.

        @param nested_list: a list that may or may not contain nested lists. Nested lists may contain additional nested
        lists.
        @type nested_list: list or list of lists
        @return: deep copy of nested_list
        @rtype: list or list of lists (matches input)
        """
        copy = []
        for item in nested_list:
            if isinstance(item, list):
                copy.append(Utilities.deep_copy(item))
            else:
                copy.append(item)
        return copy



class VariantRemapper:
    """
    Class containing functions for remapping of variants from HGVS to genomic coordinate notation.
    """

    def __init__(self):
        """
        Initializes services required for variant remapping.
        """

        # Mutalyzer webservices API URL (online tool used for remapping), See https://mutalyzer.nl/webservices for
        # documentation of all available functions.
        url = 'https://mutalyzer.nl/services/?wsdl'
        client = Client(url, cache=None)
        self.mutalyzer_interface = client.service

    def remap_variant(self, variant):
        """
        Converts a single variant provided in HGVS notation to genomic coordinate notation. Not meant for large numbers
        of variants, see submit_variant_batch for this purpose.
        See U(https://humgenprojects.lumc.nl/trac/mutalyzer/wiki/PositionConverter) for more information on acceptable
        inputs and outputs, and remapping functionality.

        @param variant: HGVS description of variant, such as NM_001100.3:c.137T>C. The portion prior to the colon is
        the refseqID used as the reference for the variant The portion after the colon is an HGVS-style description
        of the mutation (a SNP from T to C at location 137 in the example above.
        @type variant: string
        @return: Variant in HG19 coordinate notation, such as NC_000001.10:g.229568620A>G, where the portion prior
        to the colon describes the chromosome number (1 in the example) and the portion after the colon is an HGVS-
        style description of the genomic variant (a SNP from A to G at location 229568620 in the example above). Empty
        string returned if the variant could not be remapped.
        @rtype: string
        """

        # HG19 is the most recent human genome version on mutalyzer
        genome_version = 'hg19'

        try:
            result = self.mutalyzer_interface.numberConversion(genome_version, variant)
            result = result[0][0]  # converts return value to string
        except suds.WebFault:
            # Variants with syntax errors are replaced with REMAPPING_ERROR
            result = ''

        return result

    def submit_variant_batch(self, variant_list):
        """
        Submit a list of variants for batch_remapping.
        Some mistakes in HGVS notation can cause batch processing to fail. In these cases, the best
        alternative is to use remap_variant on each individual variant.

        @param variant_list: list of HGVS descriptions of variant, such as NM_001100.3:c.137T>C. The portion prior to the colon is
        the refseqID used as the reference for the variant The portion after the colon is an HGVS-style description
        of the mutation (a SNP from T to C at location 137 in the example above.
        @type variant_list: list of strings
        @return: id number for batch. Can be used to monitor job progress and retrieve results.
        @rtype: integer
        """

        mutalyzer_input = VariantRemapper.get_mutalyzer_input(variant_list)

        try:
            id_number = self.mutalyzer_interface.submitBatchJob(mutalyzer_input, 'PositionConverter', 'hg19')
        except:
            # Variants cannot be processed in batch due to formatting errors
            id_number = -1
        return id_number

    @staticmethod
    def get_mutalyzer_input(variant_list):
        """
        Creates input required for mutalyzer webAPI from list of variants.

        @param variant_list: list of variants
        @type variant_list: list of strings
        @return: string in format required for mutalyzer webAPI
        @rtype: base64 encoded string
        """
        # Encode all strings in list as binary (required for base64 encoding)
        variant_list = [x.encode() for x in variant_list]

        # Make a single \n separated list for input
        mutalyzer_input = b'\n'.join(variant_list)

        # Encode as base64 and decode so input is not binary (required for mutalyzer)
        mutalyzer_input = base64.b64encode(mutalyzer_input)
        return mutalyzer_input.decode()


    def entries_remaining_in_batch(self, id_number):
        """
        Returns the number of remaining entries to be processed in batch remapping job.

        @param id_number: id number of the batch remapping job
        @type id_number: integer
        @return: number of remaining entries to be processed in batch remapping job.
        @rtype: integer
        """

        return self.mutalyzer_interface.monitorBatchJob(id_number)


    def get_batch_results(self, id_number):
        """
        Return the results from a batch remapping job.

        @param id_number: id number of the batch remapping job
        @type id_number: integer
        @return: named tuple if lists containing chromosome_number, coordinate, ref, and alt entries for each variant
        submitted in the batch job. Coordinates are remapped to HG19 coordinates.
        @rtype: namedtuple of lists
        """
        result = self.mutalyzer_interface.getBatchJob(id_number)
        result = base64.b64decode(result).decode()
        result = result.rstrip('\n')

        # Split result string into list of lists (list of data from each row)
        row_delimiter = '\n'
        column_delimiter = '\t'
        rows = result.split(row_delimiter)
        data = [column.split(column_delimiter) for column in rows]

        # Entries with no error (most) will have no entries in that column
        variant_column = Utilities.find_string_index(data[0], 'Chromosomal Variant')
        chromosomal_variants = [x[variant_column] for x in data[1:]]

        # Extract information of interest from remapping notation
        chromosome_number = [VariantRemapper.get_chromosome_number(x) for x in chromosomal_variants]
        coordinate = [VariantRemapper.get_coordinates(x) for x in chromosomal_variants]
        ref = [VariantRemapper.get_ref(x) for x in chromosomal_variants]
        alt = [VariantRemapper.get_alt(x) for x in chromosomal_variants]

        remapping_results = namedtuple('remapping_results', 'chromosome_number coordinate ref alt')
        return remapping_results(chromosome_number, coordinate, ref, alt)

    @staticmethod
    def get_chromosome_number(mapping):
        """
        Given a mapping of the form NC_######.#:g.<HGVS notation>, where NC_###### is the chromosome reference ID,
        returns tha chromosome number from the mutation mapping. For example, calling on NC_000001.2:g.524264A>C
        would return 1. The mapping is assumed to be in valid HGVS notation.

        @param mapping: string with the mapping of a mutation in HGVS notation as shown above.
        @type mapping: string
        @return: chromosome number of the mutation. Empty string if none found.
        @rtype: string
        """

        match = re.search(r'([0]+)([1-9][0-9]?)([.])', mapping)
        if match is not None:
            return match.group(2)
        else:
            return ''

    @staticmethod
    def get_coordinates(mapping):
        """
        Given a mapping of the form NC_######.#:g.<HGVS notation>, where NC_###### is the chromosome reference ID,
        returns tha coordinate from the mutation mapping. For example, calling on NC_0000001.2:g.524264A>C
        would return 524264. The mapping is assumed to be in valid HGVS notation.

        @param mapping: string with the mapping of a mutation in HGVS notation as shown above.
        @type mapping: string
        @return: coordinates of the mutation. Empty string if none found.
        @rtype: string
        """

        match = re.search(r'([g][.])([0-9]+[_]?[0-9]*)', mapping)
        if match is not None:
            return match.group(2)
        else:
            return ''

    @staticmethod
    def get_ref(mapping):
        """
        Given a mapping of the form NC_######.#:g.<HGVS notation>, where NC_###### is the chromosome reference ID,
        returns tha REF (reference) base from the mutation mapping. For example, calling on NC_0000001:g.524264A>C
        would return A. The mapping is assumed to be in valid HGVS notation.

        @param mapping: string with the mapping of a mutation in HGVS notation as shown above.
        @type mapping: string
        @return: reference base of the mutation (base before the mutation).  Empty string if none found.
        @rtype: string
        """

        match = re.search(r'([A-Z])([>])([A-Z])', mapping)
        if match is not None:
            return match.group(1)
        else:
            return ''

    @staticmethod
    def get_alt(mapping):
        """
        Given a mapping of the form NC_######.#:g.<HGVS notation>, where NC_###### is the chromosome reference ID,
        returns tha ALT (reference) base from the mutation mapping. For example, calling on NC_0000001:g.524264A>C
        would return C. The mapping is assumed to be in valid HGVS notation.

        @param mapping: string with the mapping of a mutation in HGVS notation as shown above.
        @type mapping: string
        @return: alternate base of the mutation (base after the mutation).  Empty string if none found.
        @rtype: string
        """

        match = re.search(r'([A-Z])([>])([A-Z])', mapping)
        if match is not None:
            return match.group(3)
        else:
            return ''


class LeidenDatabase:
    """
    Should not construct directly. Use make_leiden_database factory method to obtain instances of LeidenDatabase objects.

    Class providing functions to extract information about a variants listed under a specified gene on a specified LOVD
    Leiden Database installation. For example, U(http://www.dmd.nl/nmdb2/), is a particular installation for
    variants in genes associated with Muscular Dystrophy. A list of all known installations of LOVD databases can be
    found at U(http://www.lovd.nl/2.0/index_list.php).
    """

    def __init__(self, leiden_url):
        """
        Initializes a LeidenDatabase object for the specified Leiden Database URL.

        @param leiden_url: the base URL of the particular Leiden database to be used. For example, the Leiden muscular
        dystrophy pages homepage is http://www.dmd.nl/nmdb2/. This must be a valid URL to base page of database.
        """

        self.version_number = ''
        self.leiden_home_url = ''
        self.gene_id = ''
        self.ref_seq_id = ''
        self.variant_database_url = ''
        self.gene_homepage_url = ''
        self.database_soup = ''
        self.gene_homepage_soup = ''

    @staticmethod
    def extract_LOVD_version_number(leiden_url):
        """
        Extract the version number of the LOVD installation at the specified URL.

        @param leiden_url: the base URL of the particular Leiden database to be used. For example, the Leiden muscular
        dystrophy pages homepage is http://www.dmd.nl/nmdb2/. This must be a valid URL to base page of database.
        @return: LOVD version number.
        @rtype: float
        """

        html = Utilities.get_page_html(leiden_url)

        # Extract the version number from HTML
        regex = re.compile('LOVD v\.([23])\.\d')
        results = regex.search(html)
        version_number = results.group(1)
        return float(version_number)

    def get_version_number(self):
        """
        Return version number of LOVD in use for database.

        @return: version number of LOVD in use for database
        @rtype: float
        """

        return self.version_number

    def set_gene_id(self, gene_id):
        """
        Must be called before using any non-static functions. Initializes all necessary data sources for data extraction,
        resulting in long wait times for data to download over network.

        @param gene_id: a string with the Gene ID of the gene of interest. For example, ACTA1 is the gene ID for
        actin, as specified on the Leiden Muscular Dystrophy Pages at U(http://www.dmd.nl/nmdb2/home.php?)
        @raise: ValueError if gene does not exist in the specified Leiden Database
        """

        self.gene_id = gene_id

        if gene_id in self.get_available_genes():
            # Set relevant URLs for specified gene_id
            self.variant_database_url = self.get_variant_database_url()
            self.gene_homepage_url = self.get_gene_homepage_url()

            # Extract HTML and create BeautifulSoup objects for gene_id pages
            html = Utilities.get_page_html(self.variant_database_url)
            self.database_soup = BeautifulSoup(html)

            html = Utilities.get_page_html(self.gene_homepage_url)
            self.gene_homepage_soup = BeautifulSoup(html)

            # Extract RefSeq ID for gene_id's reference transcript
            self.ref_seq_id = self.get_transcript_refseqid()
        else:
            raise ValueError('Specified gene not available in Leiden Database.')

    def get_variant_database_url(self):
        """
        Must call set_gene_id prior to use.
        Constructs URL linking to the table of variant entries for the specified gene_id on the Leiden Database site.

        @return: URL linking to the table of variant entries for the specified gene_id on the
        Leiden Database site.
        @rtype: string
        """

        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_gene_homepage_url(self):
        """
        Must call set_gene_id prior to use.
        Constructs the URL linking to the homepage for the specified gene on the Leiden Database site.

        @return: URL linking to the homepage for the specified gene on the Leiden Database site.
        @rtype: string
        """

        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_available_genes(self):
        """
        Returns a list of all genes available in the Leiden Database.

        @rtype: list of strings
        @return: list of all genes available in the LeidenDatabase.
        @rtype: string
        """

        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_link(self, link_result_set):
        """
        Given a BeautifulSoup ResultSet object containing only link tags, return relevant information for the
        given link type:
        1. PUBMED links are converted to a PMID string
        2. OMIM URLs are converted to string containing the gene ID and entry number (102610#0003 for entry 0003 in
        ACTA1 (102610), for example).
        3. Other links are returned in the format [link_string]=link_url, such as [myurl]=http://www.myurl.com for a
        link to http://www.myurl.com with the link text myurl.
        4. Links with invalid HTML markup are replaced with INVALID_LINK_MARKUP

        @param link_result_set: a collection of links. All tags must be in the following format:
        <a href = "link_url"> linkText <\a> as a BeautifulSoup ResultSet object. This is the type of object returned
        by methods such as find_all in BeautifulSoup, which is a list of matches to the specified query. See
        U(http://www.crummy.com/software/BeautifulSoup/bs4/doc/#find-all) for more information on beautiful soup 4 and
        the find_all method.
        @type link_result_set: BeautifulSoup ResultSet
        @return: list of strings with one entry for each link in link_result_set
        @rtype: list of strings

        """

        link_delimiter = ';'
        result = []
        for links in link_result_set:
            link_url = links.get('href')

            # Only get the PUBMED ID for PUBMED links
            if 'pubmed' in link_url:
                result.append("PMID=" + Utilities.get_pmid(link_url))

            # Get OMIM ID for OMIM Links
            elif 'omim' in link_url:
                result.append("OMIM=" + Utilities.get_omimid(link_url))
            # Process HGVS notation
            elif links.string and 'c.' in links.string:
                result.append("".join([self.ref_seq_id, ':', Utilities.remove_times_reported(links.string)]))
            elif links.string:
                result.append("[" + links.string + "]=" + link_url)
            else:
                result.append("INVALID_LINK_MARKUP")
        return link_delimiter.join(result)

    def get_transcript_refseqid(self):
        """
        Must call set_gene_id prior to use.
        Returns the transcript refSeq ID (the cDNA transcript used as a coordinate reference denoted by NM_... entry on
        the gene homepage on the given gene_id). For example, the ACTA1 homepage is U(http://www.dmd.nl/nmdb2/home.php)
        and the RefSeq ID for the reference transcript is "NM_001100.3".

        @rtype: string
        @return: transcript refSeqID for the object's specified gene_id. Returns an empty string if no refSeq ID
        is found for the specified gene.
        """

        # Find all links on the gene homepage
        entries = self.gene_homepage_soup.find_all('a')
        for tags in entries:
            # NM_ is unique substring to RefSeq ID. If found, return text.
            if "NM_" in tags.get_text():
                return tags.get_text()
        return ""

    def get_table_headers(self):
        """
        Must call set_gene_id prior to use.
        Returns the column labels from the table of variants in the Leiden Database variant listing for the object's
        geneID from left to right. This is the first row in the table that contains the labels for columns in the table.

        @rtype: list of strings
        @return: column labels from the table of variants in the Leiden Database variant listing for the object's
        gene_id. Returned in left to right order as they appear on the Leiden Database. Empty list returned if no
        labels are found.
        """

        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_table_data(self):
        """
        Must call set_gene_id prior to use.
        Returns a list containing lists of strings (sub-lists). Each sub-list represents a row of the table data, where
        its elements are the entries for each column within the respective row.

        @return: table data from the Leiden Database. Each sub-list represents a row of the table data, where its
        elements are the entries for each column within the respective row. The order of the sub-lists contained in
        the list matches the order of rows within the table from top to bottom and individual entries are ordered
        from left to right as they appear on the Leiden Database.
        @rtype: lists of lists of strings
        """

        total_variant_count = self.get_total_variant_count()

        # Calculate the number of pages website will use to present data
        variants_per_page = 1000  # max allowed value
        total_pages = math.ceil(total_variant_count/variants_per_page)

        # Get table data from all pages
        table_data = []
        for page_number in range(1, total_pages + 1):
            table_data.extend(self.get_table_data_page_n(page_number))

        return table_data

    def get_table_data_page_n(self, page_number):
        """
        Gets the table data from a specified page of the table of variant entries. Each page number (positive integer)
        has 1000 variants. The requested page number must not exceed the number of pages required to display the
        total number of variants. For example, using page_number = 3 for a gene with 1500 variants is not valid, as
        these are only displayed on two pages.

        @param page_number: Page containing the desired table data.
        @type page_number: integer
        @return: table data from the specified page Leiden Database. Each sub-list represents a row of the table data,
        where its elements are the entries for each column within the respective row. The order of the sub-lists
        contained in the list matches the order of rows within the table from top to bottom and individual entries are
        ordered from left to right as they appear on the Leiden Database.
        @rtype: lists of lists of strings
        """

        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_total_variant_count(self):
        """
        Get the total number of variants in the database associated with the current gene. This is the total
        number of variant entries in table of variants, not the number of unique entries.

        @return: Number of variants listed for current gene
        @rtype: number
        @raise: ValueError if the number of entries could not be found on web page
        """

        # Search for sequences of digits that are four digits or longer in length.
        m = re.compile('(\d+)\sentries')
        results = m.search(self.database_soup.get_text())

        # Return entire matched sequence (PMID)
        if results is not None:
            return int(results.group(1))
        else:
            raise ValueError('Input URL did not contain 4+ digit number.')


class LOVD2Database(LeidenDatabase):
    """
    Should not construct directly. Use make_leiden_database factory method to obtain instances of LeidenDatabase objects.

    Provides LeidenDatabse interface specific to LOVD version 2. See LeidenDatabase super class for function documentation.
    """

    def __init__(self, leiden_url):

        # Call to the super class constructor
        LeidenDatabase.__init__(self, leiden_url)
        self.version_number = 2

        # Validate and set URL for specified Leiden Database
        if not '.php' in leiden_url.lower():
            self.leiden_home_url = leiden_url

            if not leiden_url.lower().endswith('/'):
                self.leiden_home_url = leiden_url + "/"
        else:
            # Remove everything from <page>.php on from the URL to create valid base URL
            m = re.compile('[a-z]+\.php')
            result = m.search(leiden_url.lower())

            if m is not None:
                self.leiden_home_url = leiden_url[0:result.start(0)]

    def get_variant_database_url(self):

        return "".join([self.leiden_home_url, 'variants.php?action=search_unique&select_db=', self.gene_id,
                        '&limit=1000'])

    def get_gene_homepage_url(self):

        return "".join([self.leiden_home_url, 'home.php?select_db=', self.gene_id])

    def get_available_genes(self):

        # Construct URL of page containing the drop-down to select various genes
        start_url = "".join([self.leiden_home_url, '?action=switch_db'])

        # Download and parse HTML from base URL
        html = Utilities.get_page_html(start_url)
        url_soup = BeautifulSoup(html)

        # Extract all options from the SelectGeneDB drop-down control
        options = url_soup.find(id='SelectGeneDB').find_all('option')

        # Return all options in the drop-down
        available_genes = []
        for genes in options:
            available_genes.append(genes['value'])
        return available_genes

    def get_table_headers(self):

        # Find all th tags on the table of variants (column labels)
        headers = self.database_soup.find_all('th')
        result = []
        for entries in headers:
            # Column label text
            h = entries.string

            # For all entries with a string value, add them to the results (filters out extraneous th tags)
            if h is not None:
                result.append(h.strip())
        return result

    def get_table_data_page_n(self, page_number):

        # TODO can this redundancy in LOVD2/LOVD3 be eliminated?
        if page_number is not 1:

            page_url = self.variant_database_url + '&page=' + str(page_number)

            html = Utilities.get_page_html(page_url)
            database_soup = BeautifulSoup(html)
        else:
            database_soup = self.database_soup

        # id specific to data table in HTML (must be unicode due to underscore)
        table_id = "".join([u'table', u'\u005F', u'data'])

        # Extract the HTML specific to the table data
        table = database_soup.find_all(id=table_id)[0].find_all('tr')

        # First row may contain a row of images for some reason. Filter out if present.
        if table[0].find('img') is not None:
            table = table[1:]

        row_entries = []
        for rows in table:
            # Difficult to exclude column label images, as they are included as a table row with no identifier. They
            # all contain images, while none of the data rows do. This allows column label row to be filtered out.
            entries = []
            for columns in rows.find_all('td'):
                # If there are any links in the cell, process them with get_link_info
                if columns.find('a') is not None:
                    link_string = self.get_link(columns.find_all('a'))
                    entries.append(link_string)
                else:
                    entries.append(columns.string)
            row_entries.append(entries)
        return row_entries


class LOVD3Database(LeidenDatabase):
    """
    Should not construct directly. Use make_leiden_database factory method to obtain instances of LeidenDatabase objects.

    Provides LeidenDatabse interface specific to LOVD version 3. See LeidenDatabase super class for function documentation.
    """

    def __init__(self, leiden_url):

        # Call to the super class constructor
        LeidenDatabase.__init__(self, leiden_url)
        self.version_number = 3

        if not leiden_url.lower().endswith('/'):
            leiden_url += '/'

        if leiden_url.lower().endswith('genes/'):
            # Remove trailing text from URL to get the common base URL for database
            base_url_end = leiden_url.find('genes/')
            leiden_url = leiden_url[0:base_url_end]

        self.leiden_home_url = leiden_url

    def get_variant_database_url(self):

        return "".join([self.leiden_home_url, 'variants/', self.gene_id, '?page_size=1000&page=1'])

    def get_gene_homepage_url(self):

        return "".join([self.leiden_home_url, 'genes/', self.gene_id, '?page_size=1000&page=1'])

    def get_available_genes(self):

        # Construct URL of page containing the drop-down to select various genes
        start_url = "".join([self.leiden_home_url, 'genes/', '?page_size=1000&page=1'])

        # Download and parse HTML from base URL
        html = Utilities.get_page_html(start_url)
        url_soup = BeautifulSoup(html)

        # Extract all gene entries from the database homepage
        table_class = 'data'
        options = url_soup.find_all('tr', class_=table_class)

        available_genes = []
        for genes in options:
            gene_string = genes.find_all('td')[0].find('a').string
            available_genes.append(gene_string)
        return available_genes

    def get_table_headers(self):
        """
        Returns the column labels from the table of variants in the Leiden Database variant listing for the object's
        geneID from left to right. This is the first row in the table that contains the labels for entries in the rest
        of the table such as DNA change, RNA change, etc.

        @rtype: list of strings
        @return: column labels from the table of variants in the Leiden Database variant listing for the object's
        gene_id. Returned in left to right order as they appear on the Leiden Database. Empty list returned if no
        labels are found.
        """

        # Find all th tags on the table of variants (column labels)
        headers = self.database_soup.find_all('th')
        result = []
        for entries in headers:
            # Column label text
            h = entries.text

            # For all entries with a string value, add them to the results (filters out extraneous th tags)
            if h is not None:
                result.append(h.strip())
        return result

    def get_table_data_page_n(self, page_number):
        """
        Returns a list containing lists of strings (sub-lists). Each sub-list represents a row of the table data, where
        its elements are the entries for each column within the respective row.

        @rtype: lists of lists of strings
        @return: table data from the Leiden Database. Each sub-list represents a row of the table data, where its
        elements are the entries for each column within the respective row. The order of the sub-lists contained in
        the list matches the order of rows within the table from top to bottom and individual entries are ordered
        from left to right as they appear on the Leiden Database.
        """

        # TODO can this redundancy in LOVD2/LOVD3 be eliminated?
        if page_number is not 1:

            page_url = self.variant_database_url + '&page=' + str(page_number)

            html = Utilities.get_page_html(page_url)
            database_soup = BeautifulSoup(html)
        else:
            database_soup = self.database_soup

        # id specific to data table in HTML (must be unicode due to underscore)
        table_class = u'data'

        # Extract the HTML specific to the table data
        table = database_soup.find_all('tr', class_=table_class)

        # First row may contain a row of images for some reason. Filter out if present.
        if table[0].find('img') is not None:
            table = table[1:]

        row_entries = []
        for rows in table:
            # Difficult to exclude column label images, as they are included as a table row with no identifier. They
            # all contain images, while none of the data rows do. This allows column label row to be filtered out.
            entries = []
            for columns in rows.find_all('td'):
                # If there are any links in the cell, process them with get_link_info
                if columns.find('a') is not None:
                    link_string = self.get_link(columns.find_all('a'))
                    entries.append(link_string)
                else:
                    entries.append(columns.string)
            row_entries.append(entries)
        return row_entries