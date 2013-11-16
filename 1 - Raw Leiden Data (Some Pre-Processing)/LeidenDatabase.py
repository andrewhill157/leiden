from bs4 import BeautifulSoup
import urllib
import re

class LeidenDatabase:
    """
    Class providing functions to extract information about a variants listed under a specified gene on the Leiden Database (U(http://www.dmd.nl/nmdb2/home.php?action=switch_db)).
    """

    __refSeqID = ""
    __variantDatabaseURL = ''
    __geneHomepageURL = ''
    __databaseSoup = ''
    __geneHomepageSoup = ''

    def __init__(self, geneID):
        """
        Initializes a LeidenDatabase object for the specified geneID.
        Raises an exception if the specified geneID is not available in the Leiden Database (as listed in the dropdown at U(http://www.dmd.nl/nmdb2/home.php?action=switch_db))
        @param geneID: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for actin, as specified on the Leiden Database homepage (linked above).
        """
        if geneID in LeidenDatabase.getAvailableGenes():
			self.__variantDatabaseURL = LeidenDatabase.__getVariantDatabaseURL(geneID)
			self.__geneHomepageURL = LeidenDatabase.__getGeneHomepageURL(geneID)

			html = urllib.urlopen(self.__variantDatabaseURL).read()
			self.__databaseSoup = BeautifulSoup(html)

			html = urllib.urlopen(self.__geneHomepageURL).read()
			self.__geneHomepageSoup = BeautifulSoup(html)

			self.__refSeqID = self.getTranscriptRefSeqID()
        else:
            raise Exception('Specified gene not available in Leiden Database.')

    @staticmethod
    def __getVariantDatabaseURL(geneID):
        """
        Constructs the URL linking to the table of variant entries for the specified geneID on the Leiden Database site.
        @param geneID: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for actin, as specified on the Leiden Database homepage (linked above). \
        The geneID is not checked against available genes on the Leiden Database, so generated URLs are not guaranteed to be valid.
        @rtype: string
        @returns: URL (contains http://) linking to the table of variant entries for the specified geneID on the Leiden Database site. The geneID is not checked against available genes \
        on the Leiden Database, to the URL is not guaranteed to be valid.
        """
        return "".join(['http://www.dmd.nl/nmdb2/variants.php?select_db=', geneID, '&action=search_unique&order=Variant%2FDNA%2CASC&limit=1000'])

    @staticmethod
    def __getGeneHomepageURL(geneID):
        """
        Constructs the URL linking to the homepage for the specified gene on the Leiden Database site.
        @param geneID: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for actin, as specified on the Leiden Database homepage (linked above). \
        The geneID is not checked against available genes on the Leiden Database, so generated URLs are not guaranteed to be valid.
        @rtype: string
        @returns: URL linking to the homepage for the specified gene on the Leiden Database site. The geneID is not checked against available genes on the Leiden Database, to the URL is not guaranteed \
        to be valid.
        """
        return "".join(['http://www.dmd.nl/nmdb2/home.php?select_db=', geneID])
    

    def __getLinks(self, linkHTML):
        """
        Given a list of HTML formatted link tags as strings, extract relevant information and return a string containing the:
        1. PUBMED links are converted to a PMID string
        2. OMIM links are ignored entirely
        3. Other links are returned as the full URL

        @param linkHTML: a list of link tags with links to be included in result test. All tags must be in the following format: <a href = "linkURL"> linkText <\a> \
        as a BeautifulSoup result object.
        """
        LINK_DELIMITER = ';'
        result = []
        for links in linkHTML:	
                linkURL = links.get('href')
                
                # Only get the PUBMED ID for PUBMED links
                if 'pubmed' in linkURL:
                        result.append(LeidenDatabase.__getPMID(linkURL))
                        
                # Ignore OMIM links completely
                elif 'omim' in linkURL:
                        pass
                # Process HVGS notation
                elif 'c.' in links.string:
                        result.append("".join([self.__refSeqID, ':',LeidenDatabase.__removeTimesReported(links.string)]))
                else:
                        result.append(links.string)
        return LINK_DELIMITER.join(result)


    @staticmethod
    def __getPMID(linkURL):
        """
        Given a URL to a publication listed at PUBMED, return a string containing the PUBMED ID of the publication.
        
        @param linkURL: URL to the publication on PUBMED. Assumed to be a valid link to a paper on PUBMED with a PMID in the URL.
        @rtype: string
        @returns: PUBMED ID associated with linkURL (as specified by the 8 digit ID included in PUBMED URLs).
        """
        m = re.search('([/])([0-9]+)', linkURL) 
        return m.group(2) #return only the digits (PMID)

    @staticmethod
    def __removeTimesReported(hgvsNotation):
        """
        If the variant enetry contains a (Reported N times) alongside the hgvs notation, returns a new string with this parenthetical removed.
        @rtype: string
        @returns: hgvsMutation notation with instances of (Reported N times)  removed
        """
        if '(Reported' in hgvsNotation:
                return hgvsNotation[0::hvgsNotation.find('(Reported') - 1]
        else:
                return hgvsNotation

    @staticmethod
    def getAvailableGenes():
        """
        Returns a list of all genes available in the Leiden Databases (as shown in the dropdown-box at http://www.dmd.nl/nmdb2/home.php?action=switch_db. The order
        and format of elements in the list matches the order of elements in the dropdown box. 
        @rtype: list of strings
        @returns: list of all genes available in the Leiden Database
        """
        availableGenes = []
        startURL = 'http://www.dmd.nl/nmdb2/home.php?action=switch_db'
        urlText = urllib.urlopen(startURL).read()
        urlSoup = BeautifulSoup(urlText)
        options = urlSoup.find(id = 'SelectGeneDB').find_all('option')

        for genes in options:
                availableGenes.append(genes['value'])
        return availableGenes


    def getTranscriptRefSeqID(self):
        """
        Returns the transcript refSeq ID (the cDNA transcript used a a coordinate reference denoted by NM_...) for the object's geneID.
        @rtype: string
        @returns: transcript refSeqID if found for the object's geneID. Returns an empty string if no refSeq ID is found.
        """
        entries = self.__geneHomepageSoup.find_all('a')
        for tags in entries:
        	if "NM_" in tags.get_text():
        		return tags.get_text()
        return ""

    def getTableHeaders(self):
        """
        Returns the column labels from the table of variants in the Leiden Database variant listing for the object's geneID from left to right.
        This is the first row in the table that contains the labels for entries in the rest of the table such as DNA change, RNA change, etc.
        @rtype: list of strings
        @returns: column labels from the table of variants in the Leiden Database variant listing for the object's geneID. Returned in left to right order as they \
        appear on the Leiden Database. Empty list returned if no labels are found.
        """
        headers = self.__databaseSoup.find_all('th')
        result = []
        for entries in headers:
                h = entries.string

                if h is not None:
                        result.append(h)
        return result
    

    def getTableData(self):
        """
        Returns a list containing lists of strings (sub-lists). Each sub-list represents a row of the table data, where its elements are the entries for each column
        within the respective row. The order of the sub-lists contained in the list matches the order of rows within the table from top to bottom and columns from left to right
        as they appear on the Leiden Database. 
        @rtype: lists of lists of strings
        @returns: table data from the Leiden Database. Each sub-list represents a row of the table data, where its elements are the entries for each column
        within the respective row. The order of the sub-lists contained in the list matches the order of rows within the table from top to bottom and columns from left to right
        as they appear on the Leiden Database. 
        """
        table = self.__databaseSoup.find_all(onmouseout="this.className = '';")
    
        rowEntries = []
        for rows in table:
                entries = []
                for columns in rows.find_all('td'):
                        # If there are any links in the cell
                        if columns.find('a') is not None:
                                linkString = self.__getLinks(columns.find_all('a'))
                                entries.append(linkString)
                        else:
                                entries.append(columns.string)
                rowEntries.append(entries)
        return rowEntries


    def getGeneName(self):
        """
        Returns the name of the gene name associated with the database as a string (as defined by the database HTML).
        @rtype: string
        @returns: the full gene name of the gene associated as defined by the gene_name id in the Leiden Database HTML (also included in the dropdown menu located:
        http://www.dmd.nl/nmdb2/home.php?action=switch_db)
        """ 
        return self.__databaseSoup.find(id = 'gene_name').text.strip()
                        
