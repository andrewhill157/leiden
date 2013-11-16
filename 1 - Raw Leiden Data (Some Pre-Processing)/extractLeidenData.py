from LeidenDatabase import *
import sys
import argparse
import traceback


def findStringIndex(stringList, searchString):
	"""
	Given a list of strings and a string to search for, returns the first index of the first element in the list that contains the search string.
	Note that the comparison is not sensitive to case or leading or trailing whitespace characters.

	@param stringList: list of strings
	@param searchString: a string to search for in elements of stringList
	@rtype: number
	@return: index of the first instance of the search string in an element of the list. Returns -1 if the search string is not found in the list.
	"""
	i = 0
	for entry in stringList:
		if searchString.upper().strip() in entry.upper().strip():
			return i
		else: 
			i = i + 1
	return -1

def saveGeneData(geneID):
	"""
	Given a geneID, saves two files: <geneID>.txt and <geneID>_MutalizerInput.txt. from the Leiden Database (U(http://www.dmd.nl/nmdb2/home.php?action=switch_db))

	1. <geneID>.txt contains the extracted table data containing variants specific to the specified geneID in the Leiden Database. Each variant is on \
	its own line and columns are separated by commas. Header labels are included as the first line of the file.

	2. <geneID>_MutalizerInput.txt contains only the DNA Change column of <geneID>.txt (one variant per line). This file can be directly input to the \
	mutalyzer batch position converter tool by LOVD (U(https://mutalyzer.nl/batchPositionConverter))

	@param geneID: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for actin, as specified on the Leiden Database homepage (linked above).
	"""
	# Constants for file delimeters
	FILE_EXTENSION = '.txt'
	ROW_DELIMITER = '\n'
	COLUMN_DELIMITER = ','

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


def printErrors(args, gene):
	"""
	Given the list of command line arguments passed to the script and a geneID, print error messages to the console.

	@param args: a parser.parse_args() object from the argparse library. Assumed to contain an argument called debug to indicate verbosity of error messages. \
	args.debug is true, full stack traces are printed for errors. A simple error message is printed otherwise.
	@param gene: a string with the geneID of the gene that generated the error during processing.
	"""	
	if args.debug:
		print "---> " + gene + ": ERROR - NOT PROCESSED. STACK TRACE: "
		tb = traceback.format_exc()
		print tb
	else:
		print "---> " + gene + ": ERROR - NOT PROCESSED. Use --debug option for more details."


"""
Command line tool for extracting data from the leiden database.
Execute script with -h option for details on arguments and description. 
"""
parser = argparse.ArgumentParser(description="Given a geneID, saves two files: <geneID>.txt and <geneID>_MutalizerInput.txt. from the Leiden Database (http://www.dmd.nl/nmdb2/home.php?action=switch_db). \
1. <geneID>.txt contains the extracted table data containing variants specific to the specified geneID in the Leiden Database. Each variant is on \
its own line and columns are separated by commas. Header labels are included as the first line of the file. \
2. <geneID>_MutalizerInput.txt contains only the DNA Change column of <geneID>.txt (one variant per line). This file can be directly input to the \
mutalyzer batch position converter tool by LOVD (https://mutalyzer.nl/batchPositionConverter)")
group = parser.add_mutually_exclusive_group()
group.add_argument("-d", "--debug", action="store_true", help="When errors are encountered, a full stack traceback is printed.")
group.add_argument("-g", "--availableGenes", action="store_true", help="A list of all available genes is printed.")
group.add_argument("-a", "--all", action="store_true", help="Extract data for all available genes in the Leiden Database.")
parser.add_argument("geneID", help="Gene ID or multiple geneIDs to retrieve from the Leiden Database.", nargs="*")

args = parser.parse_args()

try: 
	if args.availableGenes:
		print "---> IN PROGRESS..."
		print  "\n".join(LeidenDatabase.getAvailableGenes())

	elif not args.all:

		for gene in args.geneID:
			print "---> " + gene + ": IN PROGRESS..."
			saveGeneData(gene)
			print "---> " + gene + ": COMPLETE"
	else:
		for gene in LeidenDatabase.getAvailableGenes():
			print "---> " + gene + ": IN PROGRESS..."
			saveGeneData(gene)
			print "---> " + gene + ": COMPLETE"
except:
	printErrors(args, gene)
		

from LeidenDatabase import *
import sys
import argparse
import traceback


def findStringIndex(stringList, searchString):
	"""
	Given a list of strings and a string to search for, returns the first index of the first element in the list that contains the search string.
	Note that the comparison is not sensitive to case or leading or trailing whitespace characters.

	@param stringList: list of strings
	@param searchString: a string to search for in elements of stringList
	@rtype: number
	@return: index of the first instance of the search string in an element of the list. Returns -1 if the search string is not found in the list.
	"""
	i = 0
	for entry in stringList:
		if searchString.upper().strip() in entry.upper().strip():
			return i
		else: 
			i = i + 1
	return -1

def saveGeneData(geneID):
	"""
	Given a geneID, saves two files: <geneID>.txt and <geneID>_MutalizerInput.txt. from the Leiden Database (U(http://www.dmd.nl/nmdb2/home.php?action=switch_db))

	1. <geneID>.txt contains the extracted table data containing variants specific to the specified geneID in the Leiden Database. Each variant is on \
	its own line and columns are separated by commas. Header labels are included as the first line of the file.

	2. <geneID>_MutalizerInput.txt contains only the DNA Change column of <geneID>.txt (one variant per line). This file can be directly input to the \
	mutalyzer batch position converter tool by LOVD (U(https://mutalyzer.nl/batchPositionConverter))

	@param geneID: a string with the Gene ID of the gene to be extracted. For example, ACTA1 is the gene ID for actin, as specified on the Leiden Database homepage (linked above).
	"""
	# Constants for file delimeters
	FILE_EXTENSION = '.txt'
	ROW_DELIMITER = '\n'
	COLUMN_DELIMITER = ','

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


def printErrors(args, gene):
	"""
	Given the list of command line arguments passed to the script and a geneID, print error messages to the console.

	@param args: a parser.parse_args() object from the argparse library. Assumed to contain an argument called debug to indicate verbosity of error messages. \
	args.debug is true, full stack traces are printed for errors. A simple error message is printed otherwise.
	@param gene: a string with the geneID of the gene that generated the error during processing.
	"""	
	if args.debug:
		print "---> " + gene + ": ERROR - NOT PROCESSED. STACK TRACE: "
		tb = traceback.format_exc()
		print tb
	else:
		print "---> " + gene + ": ERROR - NOT PROCESSED. Use --debug option for more details."


"""
Command line tool for extracting data from the leiden database.
Execute script with -h option for details on arguments and description. 
"""
parser = argparse.ArgumentParser(description="Given a geneID, saves two files: <geneID>.txt and <geneID>_MutalizerInput.txt. from the Leiden Database (http://www.dmd.nl/nmdb2/home.php?action=switch_db). \
1. <geneID>.txt contains the extracted table data containing variants specific to the specified geneID in the Leiden Database. Each variant is on \
its own line and columns are separated by commas. Header labels are included as the first line of the file. \
2. <geneID>_MutalizerInput.txt contains only the DNA Change column of <geneID>.txt (one variant per line). This file can be directly input to the \
mutalyzer batch position converter tool by LOVD (https://mutalyzer.nl/batchPositionConverter)")
group = parser.add_mutually_exclusive_group()
group.add_argument("-d", "--debug", action="store_true", help="When errors are encountered, a full stack traceback is printed.")
group.add_argument("-g", "--availableGenes", action="store_true", help="A list of all available genes is printed.")
group.add_argument("-a", "--all", action="store_true", help="Extract data for all available genes in the Leiden Database.")
parser.add_argument("geneID", help="Gene ID or multiple geneIDs to retrieve from the Leiden Database.", nargs="*")

args = parser.parse_args()


# User has specified the available genes option, print a list of all available genes.
if args.availableGenes:
	print "---> IN PROGRESS..."
	print  "\n".join(LeidenDatabase.getAvailableGenes())

# User has specified the all option, so extract data from all genes available on the Leiden Database
elif args.all:
	for gene in LeidenDatabase.getAvailableGenes():
		try: 
			print "---> " + gene + ": IN PROGRESS..."
			saveGeneData(gene)
			print "---> " + gene + ": COMPLETE"
		except:
			printErrors(args, gene)

# The user has not specified all, process their arguments
else:
	# No arguments passed
	if len(args.fileNames) == 0:
		print "---> NO GENES PROCESSED: Must pass at least one file or use the --all option"

	# Process each gene ID the user has passed
	else:
		for gene in args.geneID:
			try:
				print "---> " + gene + ": IN PROGRESS..."
				saveGeneData(gene)
				print "---> " + gene + ": COMPLETE"
			except:
				printErrors(args, gene)
	
		

