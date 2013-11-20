from bs4 import BeautifulSoup
import urllib.request
import re


class LeidenDatabase:
    """
    Class providing functions to extract information about a variants listed under a specified gene on the Leiden
    Database (U(http://www.dmd.nl/nmdb2/home.php?action=switch_db)).
    """

    __refSeqID = ""
    __variantDatabaseURL = ''
    __geneHomepageURL = ''
    __databaseSoup = ''
    __geneHomepageSoup = ''


    def __init__(self, gene_id):
        """
        Initializes a Leiden object for the specified gene_id.
        Raises an exception if the specified gene_id is not available in the Leiden Database (as listed in the drop-down
        at U(http://www.dmd.nl/nmdb2/home.php?action=switch_db))
        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for
        actin, as specified on the Leiden Database homepage (linked above).
        """
        if gene_id in LeidenDatabase.get_available_genes():
            self.__variantDatabaseURL = LeidenDatabase.__get_variant_database_url(gene_id)
            self.__geneHomepageURL = LeidenDatabase.__get_gene_homepage_url(gene_id)

            with urllib.request.urlopen(self.__variantDatabaseURL) as url:
                html = url.read()
                self.__databaseSoup = BeautifulSoup(html)
            with urllib.request.urlopen(self.__geneHomepageURL) as url:
                html = url.read()
                self.__geneHomepageSoup = BeautifulSoup(html)

            self.__refSeqID = self.get_transcript_refseqid()
        else:
            raise Exception('Specified gene not available in Leiden Database.')


    @staticmethod
    def __get_variant_database_url(gene_id):
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
        return "".join(['http://www.dmd.nl/nmdb2/variants.php?select_db=', gene_id,
                        '&action=search_unique&order=Variant%2FDNA%2CASC&limit=1000'])


    @staticmethod
    def __get_gene_homepage_url(gene_id):
        """
        Constructs the URL linking to the homepage for the specified gene on the Leiden Database site.
        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for \
        actin, as specified on the Leiden Database homepage (linked above). The gene_id is not checked against \
        available genes on the Leiden Database, so generated URLs are not guaranteed to be valid.
        @rtype: string
        @return: URL linking to the homepage for the specified gene on the Leiden Database site. The gene_id is not \
        checked against available genes on the Leiden Database, to the URL is not guaranteed to be valid.
        """
        return "".join(['http://www.dmd.nl/nmdb2/home.php?select_db=', gene_id])


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


    @staticmethod
    def get_available_genes():
        """
        Returns a list of all genes available in the Leiden Databases (as shown in the drop-down box at
        http://www.dmd.nl/nmdb2/home.php?action=switch_db. The order and format of elements in the list matches the
        order of elements in the drop-down box.
        @rtype: list of strings
        @return: list of all genes available in the Leiden Database
        """
        available_genes = []
        start_url = 'http://www.dmd.nl/nmdb2/home.php?action=switch_db'
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
    def save_gene_data(gene_id):
        """
        Given a gene_id, saves two files: <gene_id>.txt and <gene_id>_MutalizerInput.txt. from the Leiden Database
        (U(http://www.dmd.nl/nmdb2/home.php?action=switch_db))

        1. <gene_id>.txt contains the extracted table data containing variants specific to the specified gene_id in the
        Leiden Database. Each variant is on its own line and columns are separated by commas. Header labels are included as
        the first line of the file.

        2. <gene_id>_MutalizerInput.txt contains only the DNA Change column of <gene_id>.txt (one variant per line). This
        file can be directly input to the mutalyzer batch position converter tool by LOVD
        (U(https://mutalyzer.nl/batchPositionConverter))

        @param gene_id: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for actin,
        as specified on the Leiden Database homepage (linked above).
        """
        # Constants for file delimiters
        file_extension = '.txt'
        row_delimiter = '\n'
        column_delimiter = ','

        database = LeidenDatabase(gene_id)
        filename = "".join([gene_id, file_extension])
        mutalizer_input_file = "".join([gene_id, "_MutalizerInput", file_extension])

        # write table data to file in Unicode encoding (some characters are no ASCII encodable)
        with open(filename, 'w') as f, open(mutalizer_input_file, 'w') as mutalizer:
            entries = database.get_table_data()
            file_lines = []
            headers = database.get_table_headers()

            file_lines.append(column_delimiter.join(headers))

            hgvs_mutation_column = find_string_index(headers, u'DNA\xa0change')

            for rows in entries:
                file_lines.append(column_delimiter.join(rows))
                mutalizer.write("".join([rows[hgvs_mutation_column], "\n"]))
            f.write(row_delimiter.join(file_lines))


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


    def get_transcript_refseqid(self):
        """
        Returns the transcript refSeq ID (the cDNA transcript used a a coordinate reference denoted by NM_...) for the
        object's geneID.
        @rtype: string
        @return: transcript refSeqID for the object's geneID. Returns an empty string if no refSeq ID is found.
        """
        entries = self.__geneHomepageSoup.find_all('a')
        for tags in entries:
            if "NM_" in tags.get_text():
                return tags.get_text()
        return ""

    def get_table_headers(self):
        """
        Returns the column labels from the table of variants in the Leiden Database variant listing for the object's
        geneID from left to right. This is the first row in the table that contains the labels for entries in the rest
        of the table such as DNA change, RNA change, etc.
        @rtype: list of strings
        @return: column labels from the table of variants in the Leiden Database variant listing for the object's \
        geneID. Returned in left to right order as they appear on the Leiden Database. Empty list returned if no \
        labels are found.
        """
        headers = self.__databaseSoup.find_all('th')
        result = []
        for entries in headers:
            h = entries.string

            if h is not None:
                result.append(h)
        return result

    def get_table_data(self):
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
        table = self.__databaseSoup.find_all(onmouseout="this.className = '';")

        row_entries = []
        for rows in table:
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

    def get_gene_name(self):
        """
        Returns the name of the gene name associated with the database as a string (as defined by the database HTML).
        @rtype: string
        @return: the full gene name of the gene associated as defined by the gene_name id in the Leiden Database HTML \
        (also included in the drop-down menu located: U(http://www.dmd.nl/nmdb2/home.php?action=switch_db))
        """
        return self.__databaseSoup.find(id='gene_name').text.strip()
