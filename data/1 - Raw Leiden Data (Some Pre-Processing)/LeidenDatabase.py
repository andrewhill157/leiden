from bs4 import BeautifulSoup
import urllib

class LeidenDatabase:
    """
    TODO add documentation
    """

    __refSeqID = ""
    __variantDatabaseURL = ''
    __geneHomepageURL = ''
    __databaseSoup = ''
    __geneHomepageSoup = ''

    def __init__(self, geneID):
        """
        TODO add documentation
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
        TODO add documentation
        """
        return "".join(['http://www.dmd.nl/nmdb2/variants.php?select_db=', geneID, '&action=search_unique&order=Variant%2FDNA%2CASC&limit=1000'])

    """
    TODO document
    """
    @staticmethod
    def __getGeneHomepageURL(geneID):
        """
        TODO add documentation
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
        
        @param linkURL: URL to the publication on PUBMED
        """
        temp = linkURL.rfind('/') # should use regex to make sure there are numbers after this
        return linkURL[temp+1:]

    @staticmethod
    def __removeTimesReported(hvgsNotation):
        """
        TODO add documentation
        """
        if '(Reported' in hvgsNotation:
                return hvgsNotation[0::hvgsNotation.find('(Reported') - 1]
        else:
                return hvgsNotation

    @staticmethod
    def __isTranscriptRefSeqID(tag):
        """
        TODO add documentation
        """
        if tag.name is 'td':
                print '1'
                if tag.has_attr('a'):
                        print '2'
                        return 'NM_' in tag.a.string
        return False

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
        TODO document
        """
		entries = self.__geneHomepageSoup.find_all('a')
		for tags in entries:
			if "NM_" in tags.get_text():
				return tags.get_text()
		return 'None found'
            
    def getTableHeaders(self):
        """
        TODO add documentation
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
                        
