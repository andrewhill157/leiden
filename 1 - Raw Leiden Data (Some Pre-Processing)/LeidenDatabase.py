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
        Initializes a LeidenDatabase object for the specified gene_id.
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