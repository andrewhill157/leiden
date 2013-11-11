from bs4 import BeautifulSoup
import urllib

class LeidenDatabase:
        __refSeqID = ""
        __variantDatabaseURL = ''
        __geneHomepageURL = ''
        __databaseSoup = ''
        __geneHomepageSoup = ''

        """
        TODO - make own exception definition for gene not found in DB. Document.
        """
        def __init__(self, geneID, refSeqID):
            if geneID in LeidenDatabase.getAvailableGenes():
                self.__variantDatabaseURL = LeidenDatabase.__getVariantDatabaseURL(geneID)
                self.__geneHomepageURL = LeidenDatabase.__getGeneHomepageURL(geneID)
                
                self.__refSeqID = refSeqID
                
                html = urllib.urlopen(self.__variantDatabaseURL).read()
                self.__databaseSoup = BeautifulSoup(html)

                html = urllib.urlopen(self.__geneHomepageURL).read()
                self.__geneHomepageSoup = BeautifulSoup(html)
            else:
                raise Exception('Specified gene not available in Leiden Database.')

        """
        TODO document
        """
        @staticmethod
        def __getVariantDatabaseURL(geneID):
            return "".join(['http://www.dmd.nl/nmdb2/variants.php?select_db=', geneID, '&action=search_unique&order=Variant%2FDNA%2CASC&limit=1000'])

        """
        TODO document
        """
        @staticmethod
        def __getGeneHomepageURL(geneID):
            return "".join(['http://www.dmd.nl/nmdb2/variants.php?select_db=', geneID])
        
        """
        Given a list of HTML formatted link tags as strings, extract relevant information and return a string containing the
        PUBMED links are converted to a PMID string
        OMIM links are ignored entirely
        Other links are returned as the full URL

        @params linkHTML - a list of link tags with links to be included in result test. All tags must be in the following format: <a href = "linkURL"> linkText <\a>
        as a BeautifulSoup result object.
        TODO explain BeautifulSoup object requirement better.
        """
        def __getLinks(self, linkHTML):
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

        """
        Given a URL to a publication listed at PUBMED, return a string containing the PUBMED ID of the publication.
        
        @params linkURL - URL to the publication on PUBMED
        """
        @staticmethod
        def __getPMID(linkURL):
                temp = linkURL.rfind('/') # should use regex to make sure there are numbers after this
                return linkURL[temp+1:]

        @staticmethod
        def __removeTimesReported(hvgsNotation):
                if '(Reported' in hvgsNotation:
                        return hvgsNotation[0::hvgsNotation.find('(Reported') - 1]
                else:
                        return hvgsNotation

        """
        TODO Document
        """
        @staticmethod
        def __isTranscriptRefSeqID(tag):
                if tag.name is 'td':
                        print '1'
                        if tag.has_attr('a'):
                                print '2'
                                return 'NM_' in tag.a.string
                return False
        
        """
        Returns a list of all genes available in the Leiden Databases (as shown in the dropdown-box at http://www.dmd.nl/nmdb2/home.php?action=switch_db. The order
        and format of elements in the list matches the order of elements in the dropdown box. 
        """
        @staticmethod
        def getAvailableGenes():
                availableGenes = []
                startURL = 'http://www.dmd.nl/nmdb2/home.php?action=switch_db'
                urlText = urllib.urlopen(startURL).read()
                urlSoup = BeautifulSoup(urlText)
                options = urlSoup.find(id = 'SelectGeneDB').find_all('option')

                for genes in options:
                        availableGenes.append(genes['value'])
                return availableGenes

        """
        TODO document
        """
        def getTranscriptRefSeqID(self):
                entries = self.__geneHomepageSoup.find_all('td')
                for tags in entries:
                        if tags.find('a') is not None:
                                if tags.a.has_attr('string'):
                                        if 'NM_' in tags.a.string:
                                                return tags.a.string
                return 'None found'
                
        """
        TODO document
        """
        def getTableHeaders(self):
                headers = self.__databaseSoup.find_all('th')
                result = []
                for entries in headers:
                        h = entries.string

                        if h is not None:
                                result.append(h)
                return result
        
        """
        Returns a list containing lists of strings (sub-lists). Each sub-list represents a row of the table data, where its elements are the entries for each column
        within the respective row. The order of the sub-lists contained in the list matches the order of rows within the table from top to bottom. 
        """
        def getTableData(self):
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

        """
        Returns the name of the gene name associated with the database as a string (as defined by the database HTML).
        """ 
        def getGeneName(self):
                return self.__databaseSoup.find(id = 'gene_name').text.strip()
                        
