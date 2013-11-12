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

def findStringIndex(headerEntries, columnLabel):
	i = 0
	for entry in headerEntries:
		if columnLabel in entry:
			return i
		else: 
			i = i + 1


def getChromosomeNumber(mapping):
	m = re.search('([0]+)([1-9][0-9])([.])', mapping)
	print m.group(2)
	return m.group(2)

def getCoordinates(mapping):
	m = re.search('([g][.])([0-9]+)([_])?([0-9]+)', mapping)
	return m.group(2)

def getRef(mapping):
	m = re.search('([A-Z])([>])([A-Z])', mapping)
	return m.group(1)

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







