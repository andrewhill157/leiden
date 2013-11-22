from bs4 import BeautifulSoup
import urllib.request
import re
import traceback


class LeidenDatabase:
    """
    Class providing functions to extract information about a variants listed under a specified gene on a specified LOVD2
    Leiden Database installation. For example, U(http://www.dmd.nl/nmdb2/home.php), is a particular installation for
    variants in genes associated with Muscular Dystrophy. A list of all known installations of LOVD databases can be
    found at U(http://www.lovd.nl/2.0/index_list.php).
    """
    __leidenHomeURL = ''
    __geneID = ''
    __refSeqID = ''
    __variantDatabaseURL = ''
    __geneHomepageURL = ''
    __databaseSoup = ''
    __geneHomepageSoup = ''

    def __init__(self, leiden_url):
        """
        Initializes a Leiden object for the specified Leiden Database URL.
        @param leiden_url: the base URL of the particular Leiden database to be used. For example, the Leiden muscular \
        dystrophy pages homepage is http://www.dmd.nl/nmdb2/. This must be a valid URL to base page of database. Links \
        to specific php pages can also be passed, everything from <page>.php on in the URL will be ignored.
        """
        # Validate and set URL for specified Leiden Database
        if not '.php' in leiden_url.lower():
            self.__leidenHomeURL = leiden_url

            if not leiden_url.lower().endswith('/'):
                self.__leidenHomeURL = leiden_url + "/"
        else:
            # Remove everything from <page>.php on from the URL to create valid base URL
            m = re.compile('[a-z]+\.php')
            result = m.search(leiden_url.lower())

            if m is not None:
                self.__leidenHomeURL = leiden_url[0:result.start(0)]

    def __set_gene_id(self, gene_id):
        """
        Sets which gene_id the object is associated with.
        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for
        actin, as specified on the Leiden Database homepage (linked above).
        @raise: exception if gene does not exist in the specified Leiden Database (if the name is not an option in the
        gene selection drop-down, such as the one included at U(http://www.dmd.nl/nmdb2/home.php?)
        """

        self.__geneID = gene_id

        # Check that the gene is present in the database
        if gene_id in self.get_available_genes():
            self.__variantDatabaseURL = self.__get_variant_database_url(gene_id)
            self.__geneHomepageURL = self.__get_gene_homepage_url(gene_id)

            with urllib.request.urlopen(self.__variantDatabaseURL) as url:
                html = url.read()
                self.__databaseSoup = BeautifulSoup(html)
            with urllib.request.urlopen(self.__geneHomepageURL) as url:
                html = url.read()
                self.__geneHomepageSoup = BeautifulSoup(html)

            self.__refSeqID = self.get_transcript_refseqid(gene_id)
        else:
            raise Exception('Specified gene not available in Leiden Database.')

    def __get_variant_database_url(self, gene_id):
        """
        Constructs URL linking to the table of variant entries for the specified gene_id on the Leiden Database site.
        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for
        actin, as specified on the Leiden Database homepage (linked above). The gene_id is not checked against available
        genes on the Leiden Database, so generated URLs are not guaranteed to be valid.
        @rtype: string
        @return: URL (contains http://) linking to the table of variant entries for the specified gene_id on the \
        Leiden Database site. The gene_id is not checked against available genes on the Leiden Database, to the URL \
        is not guaranteed to be valid.
        """
        return "".join([self.__leidenHomeURL, 'variants.php?action=search_unique&select_db=', gene_id,
                        '&limit=1000'])

    def __get_gene_homepage_url(self, gene_id):
        """
        Constructs the URL linking to the homepage for the specified gene on the Leiden Database site.
        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for \
        actin, as specified on the Leiden Database homepage (linked above). The gene_id is not checked against \
        available genes on the Leiden Database, so generated URLs are not guaranteed to be valid.
        @rtype: string
        @return: URL linking to the homepage for the specified gene on the Leiden Database site. The gene_id is not \
        checked against available genes on the Leiden Database, to the URL is not guaranteed to be valid.
        """
        return "".join([self.__leidenHomeURL, 'home.php?select_db=', gene_id])

    @staticmethod
    def __get_pmid(link_url):
        """
        Given a URL to a publication listed at PUBMED, return a string containing the PUBMED ID of the publication.
        
        @param link_url: URL to the publication on PUBMED. Assumed to be a valid link to a paper on PUBMED with a PMID \
        in the URL.
        @rtype: string
        @return: PUBMED ID associated with link_url (as specified by the 8 digit ID included in PUBMED URLs).
        """
        m = re.compile('\d{8}')
        results = m.search(link_url)
        return results.group()  # return only the digits (PMID)

    @staticmethod
    def __remove_times_reported(hgvs_notation):
        """
        If the variant entry contains a (Reported N times) alongside the hgvs notation, returns a new string with
        this parenthetical removed.
        @param hgvs_notation: DNA change entry from the Leiden Database table of variants
        @rtype: string
        @return: hgvs_notation with instances of (Reported N times)  removed
        """
        if '(Reported' in hgvs_notation:
            return hgvs_notation[0::hgvs_notation.find('(Reported') - 1]
        else:
            return hgvs_notation

    def get_available_genes(self):
        """
        Returns a list of all genes available in the Leiden Databases (as shown in the drop-down box at
        http://www.dmd.nl/nmdb2/home.php?action=switch_db. The order and format of elements in the list matches the
        order of elements in the drop-down box.
        @rtype: list of strings
        @return: list of all genes available in the Leiden Database
        """
        available_genes = []
        start_url = "".join([self.__leidenHomeURL, '?action=switch_db'])
        with urllib.request.urlopen(start_url) as url:
            url_text = url.read()
            url_soup = BeautifulSoup(url_text)
            options = url_soup.find(id='SelectGeneDB').find_all('option')

        for genes in options:
            available_genes.append(genes['value'])
        return available_genes

    @staticmethod
    def find_string_index(string_list, search_string):
        """
        Given a list of strings and a string to search for, returns the first index of the first element in the list that
        contains the search string. Note that the comparison is not sensitive to case or leading or trailing whitespace
        characters.

        @param string_list: list of strings
        @param search_string: a string to search for in elements of string_list
        @rtype: number
        @return: index of the first instance of the search string in an element of the list. Returns -1 if the search string
        is not found in the list.
        """
        i = 0
        for entry in string_list:
            if search_string.upper().strip() in entry.upper().strip():
                return i
            else:
                i += 1
        return -1

    @staticmethod
    def print_errors(commandline_args, gene_name):
        """
        Given the list of command line arguments passed to the script and a geneID, print error messages to the console.

        @param commandline_args: a parser.parse_args() object from the argparse library. Assumed to contain an argument \
        called debug to indicate verbosity of error messages. args.debug is true, full stack traces are printed for errors.\
        A simple error message is printed otherwise.
        @param gene_name: a string with the geneID of the gene that generated the error during processing.
        """
        if commandline_args.debug:
            print("---> " + gene_name + ": ERROR - NOT PROCESSED. STACK TRACE: ")
            tb = traceback.format_exc()
            print(tb)
        else:
            print("---> " + gene_name + ": ERROR - NOT PROCESSED. Use --debug option for more details.")

    def __get_links(self, link_html):
        """
        Given list of HTML formatted links as strings, extract relevant information and return string containing the:
        1. PUBMED links are converted to a PMID string
        2. OMIM links are ignored entirely
        3. Other links are returned as the full URL

        @param link_html: a list of link tags with links to be included in result test. All tags must be in the \
        following format: <a href = "link_url"> linkText <\a> as a BeautifulSoup result object.
        """
        link_delimiter = ';'
        result = []
        for links in link_html:
            link_url = links.get('href')

            # Only get the PUBMED ID for PUBMED links
            if 'pubmed' in link_url:
                result.append(LeidenDatabase.__get_pmid(link_url))

            # Ignore OMIM links completely
            elif 'omim' in link_url:
                pass
            # Process HGVS notation
            elif 'c.' in links.string:
                result.append("".join([self.__refSeqID, ':', LeidenDatabase.__remove_times_reported(links.string)]))
            else:
                result.append(links.string)
        return link_delimiter.join(result)

    def get_transcript_refseqid(self, gene_id):
        """
        Returns the transcript refSeq ID (the cDNA transcript used a a coordinate reference denoted by NM_...) for the
        object's geneID.
        @rtype: string
        @return: transcript refSeqID for the object's geneID. Returns an empty string if no refSeq ID is found.
        """
        # Setup the specified gene if it has not already been done (saves reloading pages every function call)
        if gene_id is not self.__geneID:
            self.__set_gene_id(gene_id)

        entries = self.__geneHomepageSoup.find_all('a')
        for tags in entries:
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
        geneID. Returned in left to right order as they appear on the Leiden Database. Empty list returned if no \
        labels are found.
        """

        # Setup the specified gene if it has not already been done (saves reloading pages every function call)
        if gene_id is not self.__geneID:
            self.__set_gene_id(gene_id)

        headers = self.__databaseSoup.find_all('th')
        result = []
        for entries in headers:
            h = entries.string

            if h is not None:
                result.append(h)
        return result

    def get_table_data(self, gene_id):
        """
        Returns a list containing lists of strings (sub-lists). Each sub-list represents a row of the table data, where
        its elements are the entries for each column within the respective row. The order of the sub-lists contained in
        the list matches the order of rows within the table from top to bottom and columns from left to right as they
        appear on the Leiden Database.
        @rtype: lists of lists of strings
        @return: table data from the Leiden Database. Each sub-list represents a row of the table data, where its \
        elements are the entries for each column within the respective row. The order of the sub-lists contained in \
        the list matches the order of rows within the table from top to bottom and columns from left to right as they \
        appear on the Leiden Database.
        """

        # Setup the specified gene if it has not already been done (saves reloading pages every function call)
        if gene_id is not self.__geneID:
            self.__set_gene_id(gene_id)

        # id specific to data table in HTML (must be unicode due to underscore)
        table_id = "".join([u'table', u'\u005F', u'data'])

        # Extract the HTML specific to the table data
        table = self.__databaseSoup.find_all(id=table_id)[0].find_all('tr')

        # First row may contain a row of images for some reason. Filter out if present.
        if table[0].find('img') is not None:
            table = table[1:]

        row_entries = []
        for rows in table:
            # Difficult to exclude column label images, as they are included as a table row with no identifier. They
            # all contain images, while none of the data rows do. This allows column label row to be filtered out.
            entries = []
            for columns in rows.find_all('td'):
                # If there are any links in the cell
                if columns.find('a') is not None:
                    link_string = self.__get_links(columns.find_all('a'))
                    entries.append(link_string)
                else:
                    entries.append(columns.string)
            row_entries.append(entries)
        return row_entries

    def get_gene_name(self, gene_id):
        """
        Returns the name of the gene name associated with the database as a string (as defined by the database HTML).
        @rtype: string
        @return: the full gene name of the gene associated as defined by the gene_name id in the Leiden Database HTML \
        (also included in the drop-down menu located: U(http://www.dmd.nl/nmdb2/home.php?action=switch_db))
        """

        # Setup the specified gene if it has not already been done (saves reloading pages every function call)
        if gene_id is not self.__geneID:
            self.__set_gene_id(gene_id)

        return self.__databaseSoup.find(id='gene_name').text.strip()
