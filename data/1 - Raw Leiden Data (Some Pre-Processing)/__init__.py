from extractData import *
import sys

# Constant for file type
FILE_EXTENSION = '.txt'
ROW_DELIMITER = '\n'
COLUMN_DELIMITER = ','

"""
Given a list of strings and a string to search for, returns the first index of the first element in the list that contains the search string.
Note that the comparison is not sensitive to case or leading or trailing whitespace characters.
@params stringList - list of strings
@params searchString - a string to search for in elements of stringList
@returns - index of the first instance of the search string in an element of the list. Returns -1 if the search string is not found.
"""
def findStringIndex(stringList, searchString):
	i = 0
	for entry in stringList:
		if searchString.upper().strip() in entry.upper().strip():
			return i
		else: 
			i = i + 1
	return -1

def saveGeneData(geneID):
	database = LeidenDatabase(geneID)
	fileName = "".join([geneID, FILE_EXTENSION])
	mutalizerInputFile = "".join([geneID, "_MutalizerInput", FILE_EXTENSION]);

	# write table data to file in Unicode encoding (some characters are no ASCII encodable)
	with open(fileName, 'w') as f, open(mutalizerInputFile, 'w') as mutalizer:
		entries = database.getTableData()
		fileLines = []
		headers = database.getTableHeaders()

		fileLines.append(COLUMN_DELIMITER.join(headers))
		
		HVGSMutationColumn = findStringIndex(headers, u'DNA\xa0change')

		for rows in entries:
			fileLines.append(COLUMN_DELIMITER.join(rows))
			mutalizer.write("".join([rows[HVGSMutationColumn].encode("UTF-8"), "\n"]))
		f.write(ROW_DELIMITER.join(fileLines).encode("UTF-8"))
		
"""
TODO document
"""
arguments = len(sys.argv)
if arguments > 2:
	raise Exception("Must use two input arguments")
elif arguments == 2:
	saveGeneData(str(sys.argv[1]))
else: 
	for gene in LeidenDatabase.getAvailableGenes():
		try:
			saveGeneData(gene)
		except:
			print ": ".join(["Error processing gene", gene])
		

