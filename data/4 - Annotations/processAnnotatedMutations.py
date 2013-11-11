import re

# Map to go from three amino acids to one
oneLetterAACodes ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}

# Map to go from one letter amino acids to three
threeLetterAACodes = dict([[v,k] for k,v in oneLetterAACodes.items()])

def mapAACodes(code):
	if "*" in code:
		return code
	codeLength = len(code)
	if codeLength is 1:
		return threeLetterAACodes[code]
	elif codeLength is 3:
		return oneLetterAACodes[code]
	else:
		print("ERROR - only one and three letter codes supported")
		return -1

def getAnnotationInfo(line):
	# Lines are tab separated and info is in 7th column
	infoColumn = line.split("\t")[7] 
	return infoColumn.split(";")

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

		
# AA_CHANGE=-,H/D, H/D; is an example of an AA_CHANGE entry. Note that values are repeated and sometimes come up as -. Only want entries with #/# format
# LAA_CHANGE=p.(Asp294Val); is what an entry might look like for LAA_CHANGE. Can use regex to get rid of number and p.() notation to get three letter codes.
# Note that all codes should be either upper-case or lower case for comparison.
# There are some special cases where a * is used in notation, no AA_CHANGE entry appears, or LAA_CHANGE is simply a -. Need to check and make sure can account for all of these.
# TODO
def getLAAChange(annotationInfo):
	LAA_CHANGE_TAG = "LAA_CHANGE="
	for entry in annotationInfo:
		if LAA_CHANGE_TAG in entry:
			rawAnnotation = removeTag(entry, LAA_CHANGE_TAG)

			if "-" in rawAnnotation or "?" in rawAnnotation or "=" in rawAnnotation:
				return "" # no change or unknown change in LAA
			else:
				rawAnnotation = removeParentheses(rawAnnotation);
				return re.split("[0-9]+", rawAnnotation)
				

def getAAChange(annotationInfo):
	LAA_CHANGE_TAG = "LAA_CHANGE="
	AA_CHANGE_TAG = "AA_CHANGE="
	for entry in annotationInfo:
		if AA_CHANGE_TAG in entry and not LAA_CHANGE_TAG in entry:
			rawAnnotation = removeTag(entry, AA_CHANGE_TAG)
			rawAnnotation = rawAnnotation.split(",")

			for change in rawAnnotation:
				if "/" in change:
					return change.split("/")
			return ""

def isConcordant(LAA_CHANGE, AA_CHANGE):
	#TODO fix cases where this is breaking... ONLY NOT WORKING FOR SMALL SUBSET OF MUTATIONS
	LAA_CHANGE = [x.upper() for x in LAA_CHANGE]
	AA_CHANGE = [x.upper() for x in AA_CHANGE]
	LAA_CHANGE = [mapAACodes(x) for x in LAA_CHANGE]
	return LAA_CHANGE == AA_CHANGE
	
"""
This script will write all entries from the original file to a new file with a set of flags added to the beginning of each line
as a comma separated list.

Tag Definitions:
NOT_CONCORDANT - the amino acid changes caused by the mutation, noted in LAA_CHANGE and AA_CHANGE, do not match.
NO_PROTEIN_CHANGE - no protein change was caused by the mutation.
NON_SYNONYMOUS - the protein change is non-synonymous. 
NOT_HGMD - the mutation does not appear in the HGMD. 
	Potential causes are a very new publication or lack of pathogenicity.
HIGH_FREQUENCY - allele frequency above 0.5% in MacArthur Lab data sets.
	For certain genes, high allele frequencies can mean that the mutations are not worth looking into as they are common.
. - Entries that do not meet criteria for a given flag, will have a "." in place of the text for the given flag
"""
annotationFile = 'acta1_FINAL.vcf'
resultsFile = "".join(["FLAGGED_", annotationFile])
with open(annotationFile) as entries, open(resultsFile, 'w') as results:
	for line in entries:
		# Only process lines after the header
		if "#" not in line:
			flags = []
			# Step 1 - Check that LAA_CHANGE= and AA_CHANGE= entries are the same. This implies concordance. 
			annotationInfo = getAnnotationInfo(line)
			LAA_CHANGE = getLAAChange(annotationInfo)
			AA_CHANGE = getAAChange(annotationInfo)
			flags.append("".join(LAA_CHANGE))
			if LAA_CHANGE is "" or AA_CHANGE is "":
				flags.append("NO_PROTEIN_CHANGE")
			elif not isConcordant(LAA_CHANGE, AA_CHANGE):
				flags.append("NOT_CONCORDANT")
			else:
				flags.append(".")
				
				# Step 1.1 - If there is no protein change, why? Monkol expects that SEVERE_IMPACT= will be NON_SYNONYMOUS
					# TODO
			
			# Step 3 - Check for overlap with HGMD via presence of HGMD_MUT entry
			if"HGMD_MUT" not in line:
				flags.append("NOT_HGMD")
			else: 
				flags.append(".")
					
			# Step 4 - Check for high allele frequencies (not for ACTA1, but for DYSF). AC/AN is the frequency and high is > 0.5%
			# There are entries for each of the population sub-sets as well as total AC and AN values across all populations.
			# Need to look for AC_EUR, AN_EUR, etc entries, extract the values and take the ratio, flag if high. 
			# TODO
			#	flags.append("HIGH_FREQUENCY")
				
			# Write out lines to new file with [LIST_OF_FLAGS] inserted at the beginning
			results.write(",".join(flags))
			results.write(line)
			results.write("\n")
			