from bs4 import BeautifulSoup
import urllib.request
import re
import traceback
from suds.client import Client

# TODO describe leiden_url better
def get_leiden_database(leiden_url):
    """
    Factory method to generate the appropriate LeidenDatabase object depending on detected version number at specified
    URL. Only LOVD2 and LOVD3 installations are supported.
    @param leiden_url: the base URL of the particular Leiden database to be used. For example, the Leiden muscular \
    dystrophy pages homepage is http://www.dmd.nl/nmdb2/. This must be a valid URL to base page of database.
    @return: Instance of LeidenDatabase class (see definition for available methods)
    @raise: Exception if the LOVD version number is something other than 2 or 3 (unsupported)
    """
    version = get_lovd_version(leiden_url)

    if version == 2:
        database = LOVD2Database(leiden_url)
    elif version == 3:
        database = LOVD3Database(leiden_url)
    else:
        raise Exception("Unrecognized version number: " + str(version) + "!")

    return database


def get_lovd_version(leiden_url):
    """
    TODO Document
    @param leiden_url:
    @return:
    """
    with urllib.request.urlopen(leiden_url) as url:
        html = url.read().decode('utf-8')

    m = re.compile('LOVD v\.([23])\.\d')
    results = m.search(html)
    version = results.group(1)
    return float(version)


class Remapping:
    """
    TODO document
    """

    def __init__(self):
        """
        TODO document
        @return:
        """
        # Mutalyzer webservices API URL
        url = 'https://mutalyzer.nl/services/?wsdl'
        client = Client(url, cache=None)
        self.mutalyzer = client.service

    def remap_variant(self, variant):
        """
        TODO document and complete
        @param variant:
        @return:
        """

        # Remap variant to most recent genome build
        genome_build = 'hg19'
        try:
            if 'c.=' not in variant:
                result = self.mutalyzer.numberConversion(genome_build, variant)
                result = result[0][0]  # converts return value to a string
            else:
                result = 'REMAPPING_ERROR'
        except:
            result = 'REMAPPING_ERROR'

        return result


class LeidenDatabase:
    """
    Class providing functions to extract information about a variants listed under a specified gene on a specified LOVD
    Leiden Database installation. For example, U(http://www.dmd.nl/nmdb2/home.php), is a particular installation for
    variants in genes associated with Muscular Dystrophy. A list of all known installations of LOVD databases can be
    found at U(http://www.lovd.nl/2.0/index_list.php).
    """

    def __init__(self, leiden_url):
        """
        Initializes a LeidenDatabase object for the specified Leiden Database URL.

        @param leiden_url: the base URL of the particular Leiden database to be used. For example, the Leiden muscular \
        dystrophy pages homepage is http://www.dmd.nl/nmdb2/. This must be a valid URL to base page of database. Links \
        to specific php pages can also be passed, everything from <page>.php on in the URL will be ignored. \
        """
        self.version_number = ''
        self.leiden_home_url = ''
        self.gene_id = ''
        self.ref_seq_id = ''
        self.variant_database_url = ''
        self.gene_homepage_url = ''
        self.database_soup = ''
        self.gene_homepage_soup = ''

    def get_version_number(self):
        """
        TODO document
        @return:
        """
        return self.version_number

    def set_gene_id(self, gene_id):
        """
        Sets which gene_id for the object. Sets all HTML and BeautifulSoup HTML parsing objects for specified gene.

        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for
        actin, as specified on the Leiden Database homepage (linked above).
        @raise: exception if gene does not exist in the specified Leiden Database (if the name is not an option in the
        gene selection drop-down, such as the one included at U(http://www.dmd.nl/nmdb2/home.php?)
        """

        self.gene_id = gene_id

        # Check that the gene is present in the database, raise exception if not
        if gene_id in self.get_available_genes():
            # Set relevant URLs for specified gene_id
            self.variant_database_url = self.get_variant_database_url(gene_id)
            self.gene_homepage_url = self.get_gene_homepage_url(gene_id)

            # Extract HTML and create BeautifulSoup objects for gene_id pages
            with urllib.request.urlopen(self.variant_database_url) as url:
                html = url.read()
                self.database_soup = BeautifulSoup(html)
            with urllib.request.urlopen(self.gene_homepage_url) as url:
                html = url.read()
                self.gene_homepage_soup = BeautifulSoup(html)

            # Extract RefSeq ID for gene_id's reference transcript
            self.ref_seq_id = self.get_transcript_refseqid(gene_id)
        else:
            raise Exception('Specified gene not available in Leiden Database.')

    def get_variant_database_url(self, gene_id):
        """
        Constructs URL linking to the table of variant entries for the specified gene_id on the Leiden Database site.
        Note that this only constructs a valid URL for LOVD2 installations at this time.

        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for
        actin, as specified on the Leiden Muscular Dystrophy pages. The gene_id gene_id is required to match the \
        gene_id on the website exactly to generate a valid URL.
        @rtype: string
        @return: URL (contains http://) linking to the table of variant entries for the specified gene_id on the \
        Leiden Database site. Not guaranteed to be valid if the gene_id does not match the gene_id on the Leiden \
        Database site exactly.
        """
        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_gene_homepage_url(self, gene_id):
        """
        Constructs the URL linking to the homepage for the specified gene on the Leiden Database site.
        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for \
        actin, as specified on the Leiden Database homepage (linked above). The gene_id is not checked against \
        available genes on the Leiden Database, so generated URLs are not guaranteed to be valid.
        @rtype: string
        @return: URL linking to the homepage for the specified gene on the Leiden Database site. The gene_id is not \
        checked against available genes on the Leiden Database, to the URL is not guaranteed to be valid.
        """
        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    @staticmethod
    def get_pmid(link_url):
        """
        Given a URL to a publication listed on PUBMED, return a string containing the PUBMED ID of the publication.
        
        @param link_url: URL to the publication on PUBMED. Assumed to be a valid link to a paper on PUBMED. \
        For example, U(http://www.ncbi.nlm.nih.gov/pubmed/19562689) is a valid pubmed publication URL. The url must \
        contain the PMID in the URL (19562689 in the example here).
        @rtype: string
        @return: PUBMED ID associated with link_url (as specified by the N digit ID included in PUBMED URLs). Assumes \
        that PMIDs are at least 4 digits long and that no other 4 digit or longer numeric sequences are contained in \
        link_url.
        """

        # Search for sequences of digits that are four digits or longer in length.
        m = re.compile('\d{4,}')
        results = m.search(link_url)
        return results.group()  # return only the matched sequence of digits (PMID)

    @staticmethod
    def get_omimid(link_url):
        """
        Given a URL to an entry on OMIM, return a string containing the OMIM IDs for the entry.

        @param link_url: URL to the entry on OMIM. Assumed to be a valid link to an entry on PUBMED. \
        For example, U(http://www.omim.org/entry/102610#0003) is a valid link to an OMIM entry on the ACTA1 gene. \
        The url must contain the gene ID and entry number in the URL separated by a hash mark (such as, 102610#0003 in \
        the example URL.
        @rtype: string
        @return: OMIM entry associated with the URL. This consists of the gene ID (such as 102610 for ACTA1 and a \
        specific entry number (0003) separated by a hash mark (102610#0003 in the example above).
        """
        # Search for sequences of digits that are four digits or longer in length.
        m = re.compile('\d+#\d+')
        results = m.search(link_url)
        return results.group()  # return only the matched sequence of digits (PMID)

    @staticmethod
    def remove_times_reported(hgvs_notation):
        """
        If the variant entry contains (Reported N times) text alongside the HGVS variant description, returns a new \
        string with (Reported N times) removed. String is unchanged if no such substring is found.

        @param hgvs_notation: Text from entry in the DNA Change column from the Leiden Database table of variants.
        @rtype: string
        @return: hgvs_notation with instances of (Reported N times) removed. Assumes that this substring appears after \
        the HGVS notation in the entry.
        """

        if '(Reported' in hgvs_notation:
            # Return string with (Reported N times) removed
            return hgvs_notation[0::hgvs_notation.find('(Reported') - 1]
        else:
            # Return original string
            return hgvs_notation

    def get_available_genes(self):
        """
        Returns a list of all genes available in the Leiden Databases (as illustrated in the drop-down box at
        U(http://www.dmd.nl/nmdb2/home.php?action=switch_db) on the Leiden Muscular Dystrophy pages, for example). \
        The order and format of elements in the list matches the gene_ids provided in the drop-down box. See return \
        specification for more details on return format.

        @rtype: list of strings
        @return: list of all genes available in the Leiden Database associated with the object. Note that the returned \
        list only contains the gene_id values for each gene and not the full gene description provided in the \
        drop-down. For example, if the entry ACTA1 (Actin, Alpha 1 (skeletal muscle)) is listed in the drop-down box, \
        ACTA1 will be the respective entry in the returned list.
        """
        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    @staticmethod
    def find_string_index(string_list, search_string):
        """
        Given a list of strings and a string to search for, returns the first index of the first element in the list
        that contains the search string. Note that the comparison is not sensitive to case or leading or trailing
        whitespace characters.

        @param string_list: list of strings
        @param search_string: a string to search for in elements of string_list
        @rtype: number
        @return: index of the first instance of the search string in an element of the list. Returns -1 if the search \
        string is not found in the list.
        """

        i = 0
        for entry in string_list:
            if search_string.upper().strip() in entry.upper().strip():
                return i
            else:
                i += 1
        return -1

    def get_errors(self, commandline_args):
        """
        Given the list of command line arguments passed to the script, return error messages with optional stack trace.

        @param commandline_args: a parser.parse_args() object from the argparse library. Assumed to contain an \
        argument called debug to indicate verbosity of error messages. args.debug is true, full stack traces are \
        printed for errors. A simple error message is printed otherwise. See return for details.
        @rtype: string
        @return: Error message in the format "<gene_id>: ERROR - NOT PROCESSED. Use --debug option for more details." \
        if the debug option is not specified and in the format \
        "<gene_id>: ERROR - NOT PROCESSED. STACK TRACE: \n <stack trace>" if the --debug option is specified.
        """

        # Check whether the debug option has been specified in command line arguments
        if commandline_args.debug:
            tb = traceback.format_exc()
            return "".join(["---> ", self.gene_id, ": ERROR - NOT PROCESSED. STACK TRACE: \n", tb])
        else:
            return "".join(["---> ", self.gene_id, ": ERROR - NOT PROCESSED. Use --debug option for more details."])

    def get_link_info(self, link_html):
        """
        Given a BeautifulSoup ResultSet object containing only link tags, return relevant information for the \
        given link type:
        1. PUBMED links are converted to a PMID string
        2. OMIM URLs are converted to string containing the gene ID and entry number (102610#0003 for entry 0003 in \
        ACTA1 (102610), for example).
        3. Other links are returned in the format [link_string]=link_url, such as [myurl]=http://www.myurl.com for a \
        link to http://www.myurl.com with the link text myurl.
        4. Links with invalid HTML markup are replaced with INVALID_LINK_MARKUP

        @param link_html: a list of link tags with links to be included in return list. All tags must be in the \
        following format: <a href = "link_url"> linkText <\a> as a BeautifulSoup ResultSet object. This is the type \
        of object returned by methods such as find_all, which is a list of matches to the specified query. See \
        U(http://www.crummy.com/software/BeautifulSoup/bs4/doc/#find-all) for more information on beautiful soup 4 and \
        the find_all method.
        @rtype: list of strings
        @return: list of strings with one entry for each link included in the HTML, string included for each link \
        adheres to the behavior listed above.
        """

        link_delimiter = ';'
        result = []
        for links in link_html:
            link_url = links.get('href')

            # Only get the PUBMED ID for PUBMED links
            if 'pubmed' in link_url:
                result.append("PMID=" + LeidenDatabase.get_pmid(link_url))

            # Get OMIM ID for OMIM Links
            elif 'omim' in link_url:
                result.append("OMIM=" + LeidenDatabase.get_omimid(link_url))
            # Process HGVS notation
            elif links.string and 'c.' in links.string:
                result.append("".join([self.ref_seq_id, ':', LeidenDatabase.remove_times_reported(links.string)]))
            elif links.string:
                result.append("[" + links.string + "]=" + link_url)
            else:
                result.append("INVALID_LINK_MARKUP")
        return link_delimiter.join(result)

    def get_transcript_refseqid(self, gene_id):
        """
        Returns the transcript refSeq ID (the cDNA transcript used as a coordinate reference denoted by NM_... entry on\
        the gene homepage on the given gene_id). For example, the ACTA1 homepage is U(http://www.dmd.nl/nmdb2/home.php)\
        and the RefSeq ID for the reference transcript is "NM_001100.3".

        @rtype: string
        @return: transcript refSeqID for the object's specified gene_id. Returns an empty string if no refSeq ID \
        is found for the specified gene.
        """

        # Set the specified gene if it has not already been done (saves reloading pages every function call)
        if gene_id is not self.gene_id:
            self.set_gene_id(gene_id)

        # Find all links on the gene homepage
        entries = self.gene_homepage_soup.find_all('a')
        for tags in entries:
            # NM_ is unique substring to RefSeq ID. If found, return text.
            if "NM_" in tags.get_text():
                return tags.get_text()
        return ""

    def get_table_headers(self, gene_id):
        """
        Returns the column labels from the table of variants in the Leiden Database variant listing for the object's
        geneID from left to right. This is the first row in the table that contains the labels for entries in the rest
        of the table such as DNA change, RNA change, etc.

        @rtype: list of strings
        @return: column labels from the table of variants in the Leiden Database variant listing for the object's \
        gene_id. Returned in left to right order as they appear on the Leiden Database. Empty list returned if no \
        labels are found.
        """
        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_table_data(self, gene_id):
        """
        Returns a list containing lists of strings (sub-lists). Each sub-list represents a row of the table data, where
        its elements are the entries for each column within the respective row.

        @rtype: lists of lists of strings
        @return: table data from the Leiden Database. Each sub-list represents a row of the table data, where its \
        elements are the entries for each column within the respective row. The order of the sub-lists contained in
        the list matches the order of rows within the table from top to bottom and individual entries are ordered \
        from left to right as they appear on the Leiden Database.
        """
        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_gene_name(self, gene_id):
        """
        Returns the full name of the gene name associated with the database as a string. This is the full drop-down
        entry for the given gene on the Leiden Database base URL. For example, the gene name for ACTA1 on
        U(http://www.dmd.nl/nmdb2/) is "ACTA1 (ACTin, Alpha 1 (skeletal muscle))".

        @rtype: string
        @return: the full gene name of the gene associated as defined by the gene_name id in the Leiden Database HTML \
        (also included in the drop-down menu located: U(http://www.dmd.nl/nmdb2/home.php?action=switch_db))
        """

        # Set the specified gene if it has not already been done (saves reloading pages every function call)
        if gene_id is not self.gene_id:
            self.set_gene_id(gene_id)

        # gene_name id is unique to this entry
        return self.database_soup.find(id='gene_name').text.strip()


class LOVD2Database(LeidenDatabase):
    """

    """

    def __init__(self, leiden_url):

        """
        Initializes a LeidenDatabase object for the specified Leiden Database URL.
        Note that only LOVD2 installations are supported at this time.

        @param leiden_url: the base URL of the particular Leiden database to be used. For example, the Leiden muscular \
        dystrophy pages homepage is http://www.dmd.nl/nmdb2/. This must be a valid URL to base page of database. Links \
        to specific php pages can also be passed, everything from <page>.php on in the URL will be ignored. \
        """
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

    def get_variant_database_url(self, gene_id):
        """
        Constructs URL linking to the table of variant entries for the specified gene_id on the Leiden Database site.
        Note that this only constructs a valid URL for LOVD2 installations at this time.

        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for
        actin, as specified on the Leiden Muscular Dystrophy pages. The gene_id gene_id is required to match the \
        gene_id on the website exactly to generate a valid URL.
        @rtype: string
        @return: URL (contains http://) linking to the table of variant entries for the specified gene_id on the \
        Leiden Database site. Not guaranteed to be valid if the gene_id does not match the gene_id on the Leiden \
        Database site exactly.
        """

        return "".join([self.leiden_home_url, 'variants.php?action=search_unique&select_db=', gene_id,
                        '&limit=1000'])

    def get_gene_homepage_url(self, gene_id):
        """
        Constructs the URL linking to the homepage for the specified gene on the Leiden Database site.
        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for \
        actin, as specified on the Leiden Database homepage (linked above). The gene_id is not checked against \
        available genes on the Leiden Database, so generated URLs are not guaranteed to be valid.
        @rtype: string
        @return: URL linking to the homepage for the specified gene on the Leiden Database site. The gene_id is not \
        checked against available genes on the Leiden Database, to the URL is not guaranteed to be valid.
        """

        return "".join([self.leiden_home_url, 'home.php?select_db=', gene_id])

    def get_available_genes(self):
        """
        Returns a list of all genes available in the Leiden Databases (as illustrated in the drop-down box at
        U(http://www.dmd.nl/nmdb2/home.php?action=switch_db) on the Leiden Muscular Dystrophy pages, for example). \
        The order and format of elements in the list matches the gene_ids provided in the drop-down box. See return \
        specification for more details on return format.

        @rtype: list of strings
        @return: list of all genes available in the Leiden Database associated with the object. Note that the returned \
        list only contains the gene_id values for each gene and not the full gene description provided in the \
        drop-down. For example, if the entry ACTA1 (Actin, Alpha 1 (skeletal muscle)) is listed in the drop-down box, \
        ACTA1 will be the respective entry in the returned list.
        """

        # Construct URL of page containing the drop-down to select various genes
        start_url = "".join([self.leiden_home_url, '?action=switch_db'])

        # Download and parse HTML from base URL
        with urllib.request.urlopen(start_url) as url:
            url_text = url.read()
            url_soup = BeautifulSoup(url_text)

            # Extract all options from the SelectGeneDB drop-down control
            options = url_soup.find(id='SelectGeneDB').find_all('option')

        # Return all options in the drop-down
        available_genes = []
        for genes in options:
            available_genes.append(genes['value'])
        return available_genes

    def get_table_headers(self, gene_id):
        """
        Returns the column labels from the table of variants in the Leiden Database variant listing for the object's
        geneID from left to right. This is the first row in the table that contains the labels for entries in the rest
        of the table such as DNA change, RNA change, etc.

        @rtype: list of strings
        @return: column labels from the table of variants in the Leiden Database variant listing for the object's \
        gene_id. Returned in left to right order as they appear on the Leiden Database. Empty list returned if no \
        labels are found.
        """

        # Set the specified gene if it has not already been done (saves reloading pages every function call)
        if gene_id is not self.gene_id:
            self.set_gene_id(gene_id)

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

    def get_table_data(self, gene_id):
        """
        Returns a list containing lists of strings (sub-lists). Each sub-list represents a row of the table data, where
        its elements are the entries for each column within the respective row.

        @rtype: lists of lists of strings
        @return: table data from the Leiden Database. Each sub-list represents a row of the table data, where its \
        elements are the entries for each column within the respective row. The order of the sub-lists contained in
        the list matches the order of rows within the table from top to bottom and individual entries are ordered \
        from left to right as they appear on the Leiden Database.
        """

        # Set the specified gene if it has not already been done (saves reloading pages every function call)
        if gene_id is not self.gene_id:
            self.set_gene_id(gene_id)

        # id specific to data table in HTML (must be unicode due to underscore)
        table_id = "".join([u'table', u'\u005F', u'data'])

        # Extract the HTML specific to the table data
        table = self.database_soup.find_all(id=table_id)[0].find_all('tr')

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
                    link_string = self.get_link_info(columns.find_all('a'))
                    entries.append(link_string)
                else:
                    entries.append(columns.string)
            row_entries.append(entries)
        return row_entries


class LOVD3Database(LeidenDatabase):
    # TODO document
    """

    """

    # TODO document
    def __init__(self, leiden_url):

        """
        Initializes an object for the specified LOVD3 URL.

        @param leiden_url: the base URL of the particular Leiden database to be used. For example, the Leiden muscular \
        dystrophy pages homepage is http://www.dmd.nl/nmdb2/. This must be a valid URL to base page of LOVD3 database.
        Links to specific php pages can also be passed, everything from <page>.php on in the URL will be ignored. \
        """
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

    # TODO document
    def get_variant_database_url(self, gene_id):
        """
        Constructs URL linking to the table of variant entries for the specified gene_id on the Leiden Database site.
        Note that this only constructs a valid URL for LOVD2 installations at this time.

        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for
        actin, as specified on the Leiden Muscular Dystrophy pages. The gene_id gene_id is required to match the \
        gene_id on the website exactly to generate a valid URL.
        @rtype: string
        @return: URL (contains http://) linking to the table of variant entries for the specified gene_id on the \
        Leiden Database site. Not guaranteed to be valid if the gene_id does not match the gene_id on the Leiden \
        Database site exactly.
        """

        return "".join([self.leiden_home_url, 'variants/', gene_id, '?page_size=1000&page=1'])

    # TODO document
    def get_gene_homepage_url(self, gene_id):
        """
        Constructs the URL linking to the homepage for the specified gene on the Leiden Database site.
        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for \
        actin, as specified on the Leiden Database homepage (linked above). The gene_id is not checked against \
        available genes on the Leiden Database, so generated URLs are not guaranteed to be valid.
        @rtype: string
        @return: URL linking to the homepage for the specified gene on the Leiden Database site. The gene_id is not \
        checked against available genes on the Leiden Database, to the URL is not guaranteed to be valid.
        """

        return "".join([self.leiden_home_url, 'genes/', gene_id, '?page_size=1000&page=1'])

    # TODO document
    def get_available_genes(self):
        """
        Returns a list of all genes available in the Leiden Databases (as illustrated in the drop-down box at
        U(http://www.dmd.nl/nmdb2/home.php?action=switch_db) on the Leiden Muscular Dystrophy pages, for example). \
        The order and format of elements in the list matches the gene_ids provided in the drop-down box. See return \
        specification for more details on return format.

        @rtype: list of strings
        @return: list of all genes available in the Leiden Database associated with the object. Note that the returned \
        list only contains the gene_id values for each gene and not the full gene description provided in the \
        drop-down. For example, if the entry ACTA1 (Actin, Alpha 1 (skeletal muscle)) is listed in the drop-down box, \
        ACTA1 will be the respective entry in the returned list.
        """

        # Construct URL of page containing the drop-down to select various genes
        start_url = "".join([self.leiden_home_url, 'genes/', '?page_size=1000&page=1'])

        # Download and parse HTML from base URL
        with urllib.request.urlopen(start_url) as url:
            url_text = url.read()
            url_soup = BeautifulSoup(url_text)

            # Extract all gene entries from the database homepage
            table_class = 'data'
            options = url_soup.find_all('tr', class_=table_class)

        available_genes = []
        for genes in options:
            gene_string = genes.find_all('td')[0].find('a').string
            available_genes.append(gene_string)
        return available_genes

    # TODO document
    def get_table_headers(self, gene_id):
        """
        Returns the column labels from the table of variants in the Leiden Database variant listing for the object's
        geneID from left to right. This is the first row in the table that contains the labels for entries in the rest
        of the table such as DNA change, RNA change, etc.

        @rtype: list of strings
        @return: column labels from the table of variants in the Leiden Database variant listing for the object's \
        gene_id. Returned in left to right order as they appear on the Leiden Database. Empty list returned if no \
        labels are found.
        """

        # Set the specified gene if it has not already been done (saves reloading pages every function call)
        if gene_id is not self.gene_id:
            self.set_gene_id(gene_id)

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

    # TODO document
    def get_table_data(self, gene_id):
        """
        Returns a list containing lists of strings (sub-lists). Each sub-list represents a row of the table data, where
        its elements are the entries for each column within the respective row.

        @rtype: lists of lists of strings
        @return: table data from the Leiden Database. Each sub-list represents a row of the table data, where its \
        elements are the entries for each column within the respective row. The order of the sub-lists contained in
        the list matches the order of rows within the table from top to bottom and individual entries are ordered \
        from left to right as they appear on the Leiden Database.
        """

        # Set the specified gene if it has not already been done (saves reloading pages every function call)
        if gene_id is not self.gene_id:
            self.set_gene_id(gene_id)

        # id specific to data table in HTML (must be unicode due to underscore)
        table_class = u'data'

        # Extract the HTML specific to the table data
        table = self.database_soup.find_all('tr', class_=table_class)

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
                    link_string = self.get_link_info(columns.find_all('a'))
                    entries.append(link_string)
                else:
                    entries.append(columns.string)
            row_entries.append(entries)
        return row_entries