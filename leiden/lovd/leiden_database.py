import re
import math
from bs4 import BeautifulSoup
from ..input_output import web_io
from ..lovd import utilities


def make_leiden_database(leiden_url):
    """
    Factory method that returns appropriate LeidenDatabase object for lovd version installed at specified URL.
    Only LOVD2 and LOVD3 installations are supported.

    Args:
        leiden_url(str): The base URL of the particular Leiden lovd to be used.
            For example, the Leiden muscular dystrophy pages LOVD2 homepage is http://www.dmd.nl/nmdb2/. This must be a
            valid URL to base page of lovd. For LOVD3 installations, such as the Genetic Eye Disorder (GEI) Variation
            Database, the base url will be similar to http://mseqdr.lumc.edu/GEDI/. Extensions of this URL, such as
            http://mseqdr.lumc.edu/GEDI/genes or http://mseqdr.lumc.edu/GEDI/variants should not be used.

    Returns:
        LeidenDatabase: LeidenDatabase object for specified URL.

    Raises:
        Exception: If the lovd version installed at specified URL is unsupported (anything other than version 2 or 3)

    """

    version = LeidenDatabase.extract_LOVD_version_number(leiden_url)

    # Generate instance of appropriate _LeidenDatabase subclass for installed version
    if version == 2:
        database = LOVD2Database(leiden_url)
    elif version == 3:
        database = LOVD3Database(leiden_url)
    else:
        raise Exception("Unrecognized LOVD version number: " + str(version) + "!")

    return database


class LeidenDatabase:
    """
    Should not construct directly. Use make_leiden_database factory method to obtain instances of LeidenDatabase objects.

    Class providing functions to extract information about a variants listed under a specified gene on a specified lovd
    Leiden Database installation. For example, http://www.dmd.nl/nmdb2/, is a particular installation for
    variants in genes associated with Muscular Dystrophy. A list of all known installations of lovd databases can be
    found at http://www.lovd.nl/2.0/index_list.php.

    Attributes:
        version_number (float): LOVD version number
        leiden_home_url (str): url to homepage for specified gene ID
        gene_id (str): the current gene ID
        ref_seq_id (str): ID of the reference transcript used for the current reference transcript
        variant_database_url: url to the table data for current gene ID
        database_soup (BeautifulSoup): BeautifulSoup object for HTML on variant_database_url
        gene_homepage_soup (BeautifulSoup): BeautifulSoup object for HTML on gene homepage

    """

    def __init__(self, leiden_url):
        """
        Initializes a LeidenDatabase object for the specified Leiden Database URL.

        Args:
            leiden_url (str): the base URL of the particular Leiden lovd to be used.
                For example, the Leiden muscular dystrophy pages homepage is http://www.dmd.nl/nmdb2/. This must be a valid URL to base page of lovd.

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
        Extract the version number of the lovd installation at the specified URL.

        Args:
            leiden_url (str): the base URL of the particular Leiden lovd to be used.
                For example, the Leiden muscular dystrophy pages homepage is http://www.dmd.nl/nmdb2/. This must be a
                valid URL to base page of lovd.

        Returns:
            float: lovd version number

        """

        html = web_io.get_page_html(leiden_url)

        # Extract the version number from HTML
        try:
            regex = re.compile('LOVD v\.([23])\.\d', re.IGNORECASE)
            results = regex.search(html)
            version_number = results.group(1)
        except Exception as e:
            raise ValueError('No version number detected at specified URL')

        return float(version_number)

    def get_version_number(self):
        """
        Return version number of lovd in use for lovd.

        Returns:
            float: version number of lovd in use for lovd

        """

        return self.version_number

    def set_gene_id(self, gene_id):
        """
        Must be called before using any non-static functions. Initializes all necessary data sources for data extraction,
        resulting in long wait times for data to download over network.

        Args:
            param gene_id (string): a string with the Gene ID of the gene of interest. For example, ACTA1 is the gene ID for
                actin, as specified on the Leiden Muscular Dystrophy Pages at U(http://www.dmd.nl/nmdb2/home.php?)

        Raises:
            ValueError: if gene does not exist in the specified Leiden Database

        """

        self.gene_id = gene_id

        if gene_id in self.get_available_genes():
            # Set relevant URLs for specified gene_id
            self.variant_database_url = self.get_variant_database_url()
            self.gene_homepage_url = self.get_gene_homepage_url()

            # Extract HTML and create BeautifulSoup objects for gene_id pages
            html = web_io.get_page_html(self.variant_database_url)
            self.database_soup = BeautifulSoup(html)

            html = web_io.get_page_html(self.gene_homepage_url)
            self.gene_homepage_soup = BeautifulSoup(html)

            # Extract RefSeq ID for gene_id's reference transcript
            self.ref_seq_id = self.get_transcript_refseqid()
        else:
            raise ValueError('Specified gene not available in Leiden Database.')

    def get_variant_database_url(self):
        """
        Must call set_gene_id prior to use.
        Constructs URL linking to the table of variant entries for the specified gene_id on the Leiden Database site.

        Returns:
            str: URL linking to the table of variant entries for the specified gene_id on the Leiden Database site.

        """

        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_gene_homepage_url(self):
        """
        Must call set_gene_id prior to use.
        Constructs the URL linking to the homepage for the specified gene on the Leiden Database site.

        Returns:
            str: URL linking to the homepage for the specified gene on the Leiden Database site.

        """

        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_available_genes(self):
        """
        Returns a list of all genes available in the Leiden Database.

        Returns:
            list of str: list of all genes available in the LeidenDatabase.

        """

        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_link_urls(self, link_result_set):
        """
        Given a BeautifulSoup ResultSet object containing only link tags, return the URLs as a list.

        Args:
            link_result_set (BeautifulSoup ResultSet): a collection of links.
                All tags must be passed as a BeautifulSoup ResultSet object. This is the type of object returned by
                methods such as find_all in BeautifulSoup. See http://www.crummy.com/software/BeautifulSoup/bs4/doc/#find-all
                for more information.

        Returns:
            list of str: list of URLs from link_result_set

        """
        link_delimiter = ','
        result = []
        for links in link_result_set:
            link_url = links.get('href')

            # Process HGVS notation
            if links.string and ('c.' in links.string or 'p.' in links.string):
                hgvs_notation = utilities.remove_times_reported(links.string)
                hgvs_notation = utilities.correct_hgvs_parentheses(hgvs_notation)
                result.append("".join([self.ref_seq_id, ':', hgvs_notation]))
            elif links.string:
                result.append(link_url)

        return link_delimiter.join(result)

    def get_transcript_refseqid(self):
        """
        Must call set_gene_id prior to use.
        Returns the transcript refSeq ID (denoted by NM_<ID> on the gene homepage on the given gene_id).
        For example, the ACTA1 homepage is http://www.dmd.nl/nmdb2/home.php and the RefSeq ID for the reference
        transcript is "NM_001100.3".

        Returns:
            string: transcript refSeqID for set gene. Returns an empty string if no refSeq ID is found.

        """

        # Find all links on the gene homepage
        entries = self.gene_homepage_soup.find_all('a')
        for tags in entries:
            # NM_ is unique substring to RefSeq ID. If found, return text.
            if 'NM_' in tags.get_text():
                return tags.get_text()
        return ""

    def get_table_headers(self):
        """
        Must call set_gene_id prior to use.
        Returns the column labels from the table of variants in the Leiden Database for the set gene.

        Returns:
            list of str: column labels from the table of variants in the Leiden Database variant listing for the object's gene_id.
                Returned in left to right order as they appear on the Leiden Database. Empty list returned if no labels are found.

        """

        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_table_data(self):
        """
        Must call set_gene_id prior to use.
        Returns the table of variants from the set gene on the lovd.

        Returns:
            list of list of str: table of variants from the set gene on the Leiden Database.
                1st dimension is rows, 2nd is columns

        """

        total_variant_count = self.get_total_variant_count()

        # Calculate the number of pages website will use to present data
        variants_per_page = 1000  # max allowed value
        total_pages = int(math.ceil(float(total_variant_count)/float(variants_per_page)))

        # Get table data from all pages
        table_data = []
        for page_number in range(1, total_pages + 1):
            table_data.extend(self.get_table_data_page_n(page_number))

        return table_data

    def get_table_data_page_n(self, page_number):
        """
        Must call set_gene_id prior to use.
        Returns the table data from a specified page of the table of variant entries. Each page number (positive integer)
        has 1000 variants at most. The requested page number should not exceed the number of pages required to display
        the total number of variants.

        Args:
            page_number (int): page number containing the desired table data.

        Returns:
            list of list of str: table of variants from the set gene on the Leiden Database from the nth page.

        """

        raise Exception("ABSTRACT METHOD NOT IMPLEMENTED")

    def get_total_variant_count(self):
        """
        Must call set_gene_id prior to use.
        Get the total number of variants in the lovd associated with the current gene. This is the total
        number of variant entries in table of variants, not the number of unique entries.

        Returns:
            int: Number of variants listed for current gene

        Raises:
            ValueError: if the number of entries could not be found on web page

        """

        # Search for sequences of digits that are four digits or longer in length.
        m = re.compile('(\d+)\sentries')
        results = m.search(self.database_soup.get_text())

        # Return entire matched sequence (PMID)
        if results is not None:
            return int(results.group(1))
        else:
            raise ValueError('No entries found at given URL')


class LOVD2Database(LeidenDatabase):
    """
    Should not construct directly. Use make_leiden_database factory method to obtain instances of LeidenDatabase objects.

    Provides LeidenDatabse interface specific to LOVF version 2. See LeidenDatabase super class for function documentation.
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
        html = web_io.get_page_html(start_url)
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
                # Normalize headers so all lower-case, no whitespace, no non-alphanumeric characters
                h = h.lower().strip()
                h = re.sub(re.compile('[^A-Za-z0-9]'), '_', h)
                result.append(h)
        return result

    def get_table_data_page_n(self, page_number):

        # TODO can this redundancy in LOVD2/LOVD3 be eliminated?
        if page_number is not 1:

            page_url = self.variant_database_url + '&page=' + str(page_number)

            html = web_io.get_page_html(page_url)
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
            entries = []
            for columns in rows.find_all('td'):
                # If there are any links in the cell, process them with get_link_info
                if columns.find('a') is not None:
                    link_string = self.get_link_urls(columns.find_all('a'))
                    link_string = re.sub(r'\s', '', link_string)  # ensure there is no whitespace
                    entries.append(link_string)
                else:
                    column_string = columns.string.strip()  # ensure there is no whitespace
                    column_string = re.sub(r'\s', ' ', column_string)
                    entries.append(column_string)
            row_entries.append(entries)
        return row_entries


class LOVD3Database(LeidenDatabase):
    """
    Should not construct directly. Use make_leiden_database factory method to obtain instances of LeidenDatabase objects.

    Provides LeidenDatabse interface specific to lovd version 3. See LeidenDatabase super class for function documentation.
    """

    def __init__(self, leiden_url):

        # Call to the super class constructor
        LeidenDatabase.__init__(self, leiden_url)
        self.version_number = 3

        if not leiden_url.lower().endswith('/'):
            leiden_url += '/'

        if leiden_url.lower().endswith('genes/'):
            # Remove trailing text from URL to get the common base URL for lovd
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
        html = web_io.get_page_html(start_url)
        url_soup = BeautifulSoup(html)

        # Extract all gene entries from the lovd homepage
        table_class = 'data'
        options = url_soup.find_all('tr', class_=table_class)

        available_genes = []
        for genes in options:
            gene_string = genes.find_all('td')[0].find('a').string
            available_genes.append(gene_string)
        return available_genes

    def get_table_headers(self):

        # Find all th tags on the table of variants (column labels)
        headers = self.database_soup.find_all('th')
        result = []
        for entries in headers:
            # Column label text
            h = entries.text

            # For all entries with a string value, add them to the results (filters out extraneous th tags)
            if h is not None:
                h = h.lower().strip()
                h = re.sub(re.compile('[^A-Za-z0-9]'), '_', h)
                result.append(h)
        return result

    def get_table_data_page_n(self, page_number):

        # TODO can this redundancy in LOVD2/LOVD3 be eliminated?
        if page_number is not 1:

            page_url = self.variant_database_url + '&page=' + str(page_number)

            html = web_io.get_page_html(page_url)
            database_soup = BeautifulSoup(html)
        else:
            database_soup = self.database_soup

        # id specific to data table in HTML (must be unicode due to underscore)
        table_class = u'data'

        # Extract the HTML specific to the table data
        table = database_soup.find_all('tr', class_=table_class)

        if table:
            pass
        else:
            table = database_soup.find_all('tr', class_='marked')  # some LOVD3 pages had a different class

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
                    link_string = self.get_link_urls(columns.find_all('a'))
                    link_string = re.sub(r'\s', '', link_string)  # ensure there is no whitespace
                    entries.append(link_string)
                else:
                    column_string = columns.string.strip()
                    column_string = re.sub(r'\s', ' ', column_string)  # ensure there is no non-space whitespace
                    entries.append(column_string)
            row_entries.append(entries)
        return row_entries