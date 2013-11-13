import sys
import re
import os

mappingColumnLabel = "Chromosomal Variant" # TODO this may not be necessary

"""
TODO document
"""
def combineLeidenMutalyzer(rawLeidenDataFile, mutalyzerOutputFile):
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

"""
Given a list of strings and a string to search for, returns the first index of the first element in the list that contains the search string.
Note that the comparison is case sensitive and the element at a given index must only contain the search element, not exactly match the search element.
@params stringList - list of strings
@params searchString - a string to search for in elements of stringList
@returns - index of the first instance of the search string in an element of the list. Returns -1 if the search string is not found.
"""
def findStringIndex(stringList, searchString):
	i = 0
	for entry in stringList:
		if searchString in entry:
			return i
		else: 
			i = i + 1
	return -1

"""
Given a mapping of the form NC_######.#:g.<HVGS notation>, where NC_###### is the chromosome reference ID, 
returns tha chromosome number from the mutation mapping. For example, calling on NC_000001.2:g.524264A>C 
would return 1. 
The mapping is assumed to be in valid HVGS notation.
@params mapping - string with the mapping of a mutation in HVGS notation as shown above.
@returns - chromosome number of the mutation as a string
"""
def getChromosomeNumber(mapping):
	m = re.search('([0]+)([1-9][0-9]?)([.])', mapping)
	return m.group(2)

"""
Given a mapping of the form NC_######.#:g.<HVGS notation>, where NC_###### is the chromosome reference ID, 
returns tha coordinate from the mutation mapping. For example, calling on NC_0000001.2:g.524264A>C 
would return 524264. 
The mapping is assumed to be in valid HVGS notation.
@params mapping - string with the mapping of a mutation in HVGS notation as shown above.
@returns - coordinates of the mutation as a string
"""
def getCoordinates(mapping):
	m = re.search('([g][.])([0-9]+[_]?[0-9]+)', mapping)
	return m.group(2)

"""
Given a mapping of the form NC_######.#:g.<HVGS notation>, where NC_###### is the chromosome reference ID, 
returns tha REF (reference) base from the mutation mapping. For example, calling on NC_0000001:g.524264A>C 
would return A. 
The mapping is assumed to be in valid HVGS notation.
@params mapping - string with the mapping of a mutation in HVGS notation as shown above.
@returns - reference base of the mutation (base before the mutation) as a string
"""
def getRef(mapping):
	m = re.search('([A-Z])([>])([A-Z])', mapping)
	return m.group(1)

"""
Given a mapping of the form NC_######.#:g.<HVGS notation>, where NC_###### is the chromosome reference ID, 
returns tha ALT (reference) base from the mutation mapping. For example, calling on NC_0000001:g.524264A>C 
would return C. 
The mapping is assumed to be in valid HVGS notation.
@params mapping - string with the mapping of a mutation in HVGS notation as shown above.
@returns - alternate base of the mutation (base after the mutation) as a string
"""
def getAlt(mapping):
	m = re.search('([A-Z])([>])([A-Z])', mapping)
	return m.group(3)

"""
TODO document
"""
def convertToVCF(combinedData):
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

"""
Writes a list of lists to a file with a specified file name. 
Elements of listOfLists (lists themselves) are each written to their own line in the file with a tab character separating each entry.
@params listOfLists - list whose elements are lists of data elements such as strings, ints, etc. 
@params fileName - File name of the file to write listOfLists to. Must be a valid file name with extension.
"""
def writeListOfListsToFile(listOfLists, fileName):
	with open(fileName, 'w') as f:
		for lists in listOfLists:
			f.write("\t".join(lists))

"""
Removes the file extension from a given file. 
@params fileName - name of the file whose extension is to be removed. fileName may not contain any other "." characters other than that that
preceeds the file extension. Files without extensions are returned unmodified.
@returns fileName with no file extension ("." character is also removed). For example an input of myFile.txt would return a new string, myFile
"""
def removeFileExtension(fileName):
	extensionIndex = fileName.find(".")
	return fileName[0:extensionIndex]

"""
TODO document
"""
def combineMappingData(rawLeidenDataFile):
	try:
		if ".txt" in files and "MutalizerOutput" not in files and "_MAPPED" not in files:
			rawLeidenDataFile = files
			
			mutalyzerOutputFile = "".join([removeFileExtension(files), "_MutalizerOutput.txt"])

			combinedData = combineLeidenMutalyzer(rawLeidenDataFile, mutalyzerOutputFile)
			VCFData = convertToVCF(combinedData)
			writeListOfListsToFile(VCFData, "_".join([removeFileExtension(rawLeidenDataFile), "_MAPPED.txt"]))
	except: 
		print ": ".join(["Error processing gene", files])

"""
This script can be called in two ways:

1. Call with no arguments. This will combine and process all text files in the directory that do not contain _MAPPED.
Note that this requires any number of text files with raw Leiden Database data along with a set of corresponding files with 
_MutalizerOutput as a post-fix to the original file name. These corresponding files are the output of the Mutalyzer remapping tool, 
and contain remapped coordinates for mutations (one per line). For example ACTA1.txt and ACTA1_MutalizerOutput are corresponding.
The product of a this call is a file with combined and processed data with the same name as the original raw Leiden Database data
file with a post-fix of _MAPPED (ATCA1_MAPPED.txt, for example).

2. Call with arguments specifying the name of the raw Leiden Database files. Any number of files can be included in the list, each 
separated by a space. For example, python compileMappings.py ACTA1.txt DYSF.txt, etc. The same requirements for corresponding file names
and output file names still apply when using this method. This is useful for specifying specific files to operate on rather than an entire
directory. 
"""	
arguments = len(sys.argv)
if arguments > 1:
	for files in sys.argv[1:]:
		combineMappingData(files)
else:
	for files in os.listdir(os.path.dirname(os.path.abspath(__file__))):
		combineMappingData(files)










