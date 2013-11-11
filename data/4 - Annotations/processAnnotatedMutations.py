import re
import sys
 
 
def mapAACodes(code):
	# Mapping from three letter amino acid codes to one letter amino acid codes 
	oneLetterAACodes ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
	'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
	'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
	'GLY':'G', 'PRO':'P', 'CYS':'C'}
	
	# Mapping from one letter amino acid codes to three letter amino acid codes 
	threeLetterAACodes = dict([[v,k] for k,v in oneLetterAACodes.items()])

	if "*" in code:
		return code
	codeLength = len(code)
	if codeLength is 1:
		return threeLetterAACodes[code]
	elif codeLength is 3:
		return oneLetterAACodes[code]
	else:
		return "ERROR"

def getAnnotationInfo(line):
	# Lines in original file are tab separated and info is in 8th column
	infoColumn = line.split("\t")[7] 
	return infoColumn.split(";")

def getTaggedEntry(annotationInfo, tag):
	for entry in annotationInfo:
		if entry.startswith(tag):
			return removeTag(entry, tag)
	return ""
			
def removeTag(infoText, tag):
	tagLength = len(tag)
	return infoText[tagLength:].strip()

def removeParentheses(annotationText):
	openingParenthesisIndex = annotationText.find("(");
	pdotNotationIndex = annotationText.find(".");
	
	if openingParenthesisIndex > 0:
		closingParenthesisIndex = annotationText.find(")")
		return annotationText[openingParenthesisIndex+1:closingParenthesisIndex]
	elif pdotNotationIndex > 0:
		return annotationText[pdotNotationIndex+1:]
	else: return annotationText

		
def getLAAChange(annotationInfo):
	LAA_CHANGE_TAG = "LAA_CHANGE="
	rawAnnotation = getTaggedEntry(annotationInfo, LAA_CHANGE_TAG)
	
	if "-" in rawAnnotation or "?" in rawAnnotation or "=" in rawAnnotation:
		return "" # no change or unknown change in LAA
	else:
		rawAnnotation = removeParentheses(rawAnnotation);
		return re.split("[0-9]+", rawAnnotation)
				

def getAAChange(annotationInfo):
	AA_CHANGE_TAG = "AA_CHANGE="
	rawAnnotation = getTaggedEntry(annotationInfo, AA_CHANGE_TAG)
	rawAnnotation = rawAnnotation.split(",")

	for change in rawAnnotation:
		if "/" in change:
			return change.split("/")
	return ""

def isConcordant(LAA_CHANGE, AA_CHANGE):
	try:
		if len(LAA_CHANGE) is not 2 or len(AA_CHANGE) is not 2:
			return "ERROR"
		LAA_CHANGE = [x.upper() for x in LAA_CHANGE]
		AA_CHANGE = [x.upper() for x in AA_CHANGE]
		LAA_CHANGE = [mapAACodes(x) for x in LAA_CHANGE]
		if LAA_CHANGE == AA_CHANGE:
			return "OK"
		else: 
			return "FAIL"
	except: 
		return "ERROR"

def getSevereImpact(annotationInfo):
	severeImpactTag = "SEVERE_IMPACT="
	return getTaggedEntry(annotationInfo, severeImpactTag)


def getAlleleFrequency(annotationInfo):
	ACTag = "AC_MAC26K="
	ANTag = "AN_MAC26K="
	
	AC = getTaggedEntry(annotationInfo, ACTag)
	AN = getTaggedEntry(annotationInfo, ANTag)
	
	if AC == "" or AN == "":
		return ""
	else:
		AC = float(AC)
		AN = float(AN)
		return (AC/AN)*100
		
"""
This script must be called in the following format:
	python processAnnotatedMutations.py myAnnotatedMutations.vcf
Note that the .vcf file parameter should be adjusted to the file name of the file to be processed.

The script will output a file in the same directory as the script called FLAGGED_<original_file_name>.vcf
For example, the command shown in the paragraph above would produce a file called FLAGGED_myAnnotatedMutations.vcf

The output file will contain all header same information as the original file, except the column will now contain a semi-colon delimited list of 
tags and corresponding values for each mutation: TAG1=VALUE;TAG2=VALUE;TAG3=VALUE... 
Note that the original columns in the file will all be shifted right one column in the output.

Tag and value definitions in resulting output file:
CONCORDANCE=
	- CONCORDANT: The amino acid changes predicted by LAA_CHANGE and AA_CHANGE in the input file match
	- NOT_CONCORDANT: The amino acid changes predicted by LAA_CHANGE and AA_CHANGE in the input file do not match
	- ERROR: There was an error processing the LAA_CHANGE and/or AA_CHANGE entry
	- NO_AA_CHANGE: No amino acid change is predicted by the annotation

SEVERE_IMPACT=
	- The value of SEVERE_IMPACT in the original file is simply copied as the value of this tag for convenience

HGMD=
	- OVERLAPS - An HGMD_MUT entry is present in the original file
	- NOT_IN_HGMD - No HGMD_MUT entry is present in the original file
	
ALLELE_FREQUENCY= 
	- LOW_FREQUENCY: Overall allele frequency,(AC_MAC26K/AN_MAC26K)*100, is less than or equal to 0.5%
	- HIGH_FREQUENCY: Overall allele frequency,(AC_MAC26K/AN_MAC26K)*100, is greater than 0.5%
	- NOT_IN_26K_DATABASE: There were no AC_MAC26K or AN_MAC26K entries for this mutation
"""
if len(sys.argv) is not 2:
	raise Exception("Must use one input argument")
	
annotationFile = str(sys.argv[1])
resultsFile = "".join(["FLAGGED_", annotationFile])
with open(annotationFile) as entries, open(resultsFile, 'w') as results:
	for line in entries:
		# Only process flags for lines that are not header information
		if "#" not in line:
			flags = []
			# STEP 1 - Check that LAA_CHANGE= and AA_CHANGE= entries are the same (that entries are concordant)
			annotationInfo = getAnnotationInfo(line)
			LAA_CHANGE = getLAAChange(annotationInfo)
			AA_CHANGE = getAAChange(annotationInfo)
			
			CONCORDANCE_FLAG = "CONCORDANCE="
			if LAA_CHANGE is "" or AA_CHANGE is "":
				flags.append("".join([CONCORDANCE_FLAG, "NO_AA_CHANGE"]))
				
			else:
				concordance = isConcordant(LAA_CHANGE, AA_CHANGE)
				if concordance == "OK":
					flags.append("".join([CONCORDANCE_FLAG, "CONCORDANT"]))
				elif concordance == "FAIL":
					flags.append("".join([CONCORDANCE_FLAG,"NOT_CONCORDANT"]))
				elif concordance == "ERROR":
					flags.append("".join([CONCORDANCE_FLAG,"ERROR"]))	
			
			# STEP 2 - Extract the SEVERE_IMPACT entry
			SEVERE_IMPACT_FLAG = "SEVERE_IMPACT="
			severeImpact = getSevereImpact(annotationInfo)
			flags.append("".join([SEVERE_IMPACT_FLAG,severeImpact]))
				
			# STEP 3 - Check for overlap with HGMD
			HGMD_FLAG = "HGMD="
			if "HGMD_MUT" not in line:
				flags.append("".join([HGMD_FLAG,"NOT_IN_HGMD"]))
			else: 
				flags.append("".join([HGMD_FLAG,"OVERLAPS"]))
					
			# Step 4 - Check for allele frequency information (above or below overall 0.5% frequency
			ALLELE_FREQUENCY_TAG = "ALLELE_FREQUENCY="
			alleleFrequency = getAlleleFrequency(annotationInfo)
			if alleleFrequency == "":
				flags.append("".join([ALLELE_FREQUENCY_TAG,"NOT_IN_26K_DATABASE"]))
			elif alleleFrequency > 0.5:
				flags.append("".join([ALLELE_FREQUENCY_TAG,"HIGH_FREQUENCY"]))
			else: 
				flags.append("".join([ALLELE_FREQUENCY_TAG,"LOW_FREQUENCY"]))
				
			# Write out lines to new file with [LIST_OF_FLAGS] inserted as the new first column
			results.write(";".join(flags))
			results.write("\t")
			results.write(line)
		
		# Copy header data unchanged from original file
		elif line.startswith("##"):
			results.write(line)
			
		# Copy column labels with new column, FLAGS, inserted as label for first column
		elif line.startswith("#"):
			results.write("".join([line[0], "FLAGS", "\t", line[1:]]))
			