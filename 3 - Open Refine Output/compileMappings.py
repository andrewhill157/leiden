import sys
import re
import os
import argparse
import traceback

def combineLeidenMutalyzer(rawLeidenDataFile, mutalyzerOutputFile):
	"""
	TODO document
	"""
	with open(rawLeidenDataFile) as leidenData, open(mutalyzerOutputFile) as mutalyzerOutput:
		errors = []
		mappings = []

		# Store all entries in from the mutalizer output
		for line in mutalyzerOutput:
			columns = line.split("\t");
			errors.append(columns[1])
			mappings.append(columns[2]);

		# Insert the mutalyzer output into the data
		combinedData = []
		i = 0
		for line in leidenData:
			temp = line.split(",")

			temp.insert(2, errors[i])
			temp.insert(3, mappings[i])
			combinedData.append(temp)
			i = i + 1

		return combinedData


def findStringIndex(stringList, searchString):
	"""
	Given a list of strings and a string to search for, returns the first index of the first element in the list that contains the search string.
	Note that the comparison is case sensitive and the element at a given index must only contain the search element, not exactly match the search element.
	@param stringList: list of strings
	@param searchString: a string to search for in elements of stringList
	@rtype: number
	@returns: index of the first instance of the search string in an element of the list. Returns -1 if the search string is not found.
	"""
	i = 0
	for entry in stringList:
		if searchString in entry:
			return i
		else: 
			i = i + 1
	return -1


def getChromosomeNumber(mapping):
	"""
	Given a mapping of the form NC_######.#:g.<HVGS notation>, where NC_###### is the chromosome reference ID, 
	returns tha chromosome number from the mutation mapping. For example, calling on NC_000001.2:g.524264A>C 
	would return 1. 
	The mapping is assumed to be in valid HVGS notation.
	@param mapping: string with the mapping of a mutation in HVGS notation as shown above.
	@rtype: string
	@returns: chromosome number of the mutation as a string
	"""
	m = re.search('([0]+)([1-9][0-9]?)([.])', mapping)
	return m.group(2)


def getCoordinates(mapping):
	"""
	Given a mapping of the form NC_######.#:g.<HVGS notation>, where NC_###### is the chromosome reference ID, 
	returns tha coordinate from the mutation mapping. For example, calling on NC_0000001.2:g.524264A>C 
	would return 524264. 
	The mapping is assumed to be in valid HVGS notation.
	@param mapping: string with the mapping of a mutation in HVGS notation as shown above.
	@rtype: string
	@returns: coordinates of the mutation as a string
	"""
	m = re.search('([g][.])([0-9]+[_]?[0-9]+)', mapping)
	return m.group(2)


def getRef(mapping):
	"""
	Given a mapping of the form NC_######.#:g.<HVGS notation>, where NC_###### is the chromosome reference ID, 
	returns tha REF (reference) base from the mutation mapping. For example, calling on NC_0000001:g.524264A>C 
	would return A. 
	The mapping is assumed to be in valid HVGS notation.
	@param mapping: string with the mapping of a mutation in HVGS notation as shown above.
	@rtype: string
	@returns: reference base of the mutation (base before the mutation) as a string
	"""
	m = re.search('([A-Z])([>])([A-Z])', mapping)
	return m.group(1)

	
def getAlt(mapping):
	"""
	Given a mapping of the form NC_######.#:g.<HVGS notation>, where NC_###### is the chromosome reference ID, 
	returns tha ALT (reference) base from the mutation mapping. For example, calling on NC_0000001:g.524264A>C 
	would return C. 
	The mapping is assumed to be in valid HVGS notation.
	@param mapping: string with the mapping of a mutation in HVGS notation as shown above.
	@rtype: string
	@returns: alternate base of the mutation (base after the mutation) as a string
	"""
	m = re.search('([A-Z])([>])([A-Z])', mapping)
	return m.group(3)


def convertToVCF(combinedData):
	"""
	TODO document
	"""
	mappingIndex = findStringIndex(combinedData[0], "Chromosomal Variant")
	VCFMapping = []
	combinedData[0].insert(mappingIndex, "ALT")
	combinedData[0].insert(mappingIndex, "REF")
	combinedData[0].insert(mappingIndex, "COORDINATE")
	combinedData[0].insert(mappingIndex, "CHROMOSOME NUMBER")
	VCFMapping.append(combinedData[0])
	for row in combinedData[1:]:
		mapping = row[mappingIndex]

		if mapping is not "":
			chromosome = getChromosomeNumber(mapping)

			if ">" in mapping:
				coordinates = getCoordinates(mapping)
				ref = getRef(mapping)
				alt = getAlt(mapping)

			else:
				# Not yet processing anything more complicated than SNPs
				coordinates = mapping
				ref = ""
				alt = ""
		else: 
			chromosome = ""
			coordinates = ""
			ref = ""
			alt = ""

		row.insert(mappingIndex, alt)
		row.insert(mappingIndex, ref)
		row.insert(mappingIndex, coordinates)
		row.insert(mappingIndex, chromosome)
		VCFMapping.append(row)

	return VCFMapping

def writeListOfListsToFile(listOfLists, fileName):
	"""
	Writes a list of lists to a file with a specified file name. 
	Elements of listOfLists (lists themselves) are each written to their own line in the file with a tab character separating each entry.
	@param listOfLists: list whose elements are lists of data elements such as strings, ints, etc. 
	@param fileName: File name of the file to write listOfLists to. Must be a valid file name with extension.
	"""
	with open(fileName, 'w') as f:
		for lists in listOfLists:
			f.write("\t".join(lists))


def removeFileExtension(fileName):
	"""
	Removes the file extension from a given file. 
	@param fileName: name of the file whose extension is to be removed. Files without extensions are returned unmodified.
	@rtype: string
	@returns: fileName with no file extension ("." character is also removed). For example an input of myFile.txt would return a new string, myFile
	"""
	return os.path.splitext(fileName)[0]


def combineMappingData(rawLeidenDataFile):
	"""
	TODO document
	"""
	# Do not want to process MutalizerOutput files of files that have already been combined (_MAPPED)
	if ".txt" in files and "MutalizerOutput" not in files and "_MAPPED" not in files:
		rawLeidenDataFile = files
		
		mutalyzerOutputFile = "".join([removeFileExtension(files), "_MutalizerOutput.txt"])

		combinedData = combineLeidenMutalyzer(rawLeidenDataFile, mutalyzerOutputFile)
		VCFData = convertToVCF(combinedData)
		writeListOfListsToFile(VCFData, "_".join([removeFileExtension(rawLeidenDataFile), "_MAPPED.txt"]))

def getFilesInCurrentDirectory():
	"""
	Returns a list of all files in the current directory (the directory the script is running in).
	@rtype: list of strings
	@returns: list of all files in the current directory (the directory the script is running in).
	"""
	return os.listdir(os.path.dirname(os.path.abspath(__file__)))

def printErrors(args, fileName):
	"""
	Given the list of command line arguments passed to the script and a fileName, print error messages to the console.

	@param args: a parser.parse_args() object from the argparse library. Assumed to contain an argument called debug to indicate verbosity of error messages. \
	args.debug is true, full stack traces are printed for errors. A simple error message is printed otherwise.
	@param fileName: a string with the fileName of the file that generated the error during processing.
	"""	
	if args.debug:
		print "---> " + fileName + ": ERROR - NOT PROCESSED. STACK TRACE: "
		tb = traceback.format_exc()
		print tb
	else:
		print "---> " + fileName + ": ERROR - NOT PROCESSED. Use --debug option for more details."


"""
Command line tool for extracting data from the leiden database.
Execute script with -h option for details on arguments and description. 
"""
parser = argparse.ArgumentParser(description="Combines output from the Mutalyzer Batch Position Converter Tool and the original extracted data from the Leiden Database entry for a given gene. \
	This produces a file in the called <geneID>_MAPPED, which contains the combined information of <geneID>.txt (original extracted Leiden Database data) and <geneID>_MutalizerOutput (contains the \
		mutalizer output from all variants from the original Leiden Database. These two files are required and must be present in the same directory to run the script and must be named in this manner. \
		 Note that the output file is also produced in this same directory. It is assumed that all variants from the original data are present in the both files.")
group = parser.add_mutually_exclusive_group()
group.add_argument("-d", "--debug", action="store_true", help="When errors are encountered, a full stack traceback is printed.")
group.add_argument("-a", "--all", action="store_true", help="All files with a .txt extension that do not contain MutalizerOutput or _MAPPED in the filename are processed.")
parser.add_argument("fileNames", help="File name or multiple file names to compile. This should be a text file with the original extracted Leiden Database data (the file ACTA1.txt when calling \
	python extractLeidenData.py ACTA1, for example)", nargs="*")

args = parser.parse_args()

# User has specified the aal option, process all files in the directory
if args.all:
	for files in getFilesInCurrentDirectory():
		try:
			combineMappingData(files)
			print "---> " + files + ": COMPLETE"
		except:
			printErrors(args, files)

# The user has not specified all, process their arguments
else:
	# No arguments passed
	if len(args.fileNames) == 0:
		print "---> NO FILES PROCESSED: Must pass at least one file or use the --all option"
	# Process each file the user passes as an argument
	else:
		for files in args.fileNames:
			try:
				combineMappingData(files)
				print "---> " + files + ": COMPLETE"
			except:
				printErrors(args, files)










