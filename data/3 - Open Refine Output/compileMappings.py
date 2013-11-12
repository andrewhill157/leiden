import sys
import re
import os

mappingColumnLabel = "Chromosomal Variant"

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


def writeListOfListsToFile(listOfLists, fileName):
	with open(fileName, 'w') as f:
		for lists in listOfLists:
			f.write("\t".join(lists))


def removeFileExtension(fileName):
	extensionIndex = fileName.find(".")
	return fileName[0:extensionIndex]



for files in os.listdir(os.path.dirname(os.path.abspath(__file__))):
	if ".txt" in files and "MutalizerOutput" not in files and "MAPPED_" not in files:
		print files
		rawLeidenDataFile = files
		
		mutalyzerOutputFile = "".join([removeFileExtension(files), "_MutalizerOutput.txt"])

		combinedData = combineLeidenMutalyzer(rawLeidenDataFile, mutalyzerOutputFile)
		VCFData = convertToVCF(combinedData)
		writeListOfListsToFile(VCFData, "_".join([removeFileExtension(rawLeidenDataFile), "_MAPPED.txt"]))







