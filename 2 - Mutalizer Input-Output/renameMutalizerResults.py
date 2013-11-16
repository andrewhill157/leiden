import glob, os
import shutil

results = []
inputs = []
for x in sorted(glob.glob("*.txt")):
	if "Results" in x:
		results.append(x)
	elif "Output" not in x:
		inputs.append(x)

if len(results) == len(inputs):
	for i in range(0, len(results)):
		inputFile = inputs[i]
		resultFile = results[i]
		geneID = inputFile[0:inputFile.index("_")]
		shutil.copyfile(resultFile, geneID + "_MutalizerOutput.txt")
else:
	print "Number of inputs and results does not match!"

		