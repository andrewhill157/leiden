from extractData import *
import sys

# Constant for file type
FILE_EXTENSION = '.txt'
ROW_DELIMITER = '\n'
COLUMN_DELIMITER = ','

# User passes the geneID and refSeqID they want to extract
if len(sys.argv) is not 3:
	raise Exception("Must use two input arguments")
	
geneID = str(sys.argv[1])
refSeqID = str(sys.argv[2])

database = LeidenDatabase(geneID, refSeqID)
fileName = "".join([geneID, FILE_EXTENSION])
mutalizerInputFile = "".join([geneID, "_MutalizerInput", FILE_EXTENSION]);

# write table data to file in Unicode encoding (some characters are no ASCII encodable)
with open(fileName, 'w') as f, open(mutalizerInputFile, 'w') as mutalizer:
	entries = database.getTableData()
	fileLines = []
	headers = database.getTableHeaders()
	fileLines.append(COLUMN_DELIMITER.join(headers))
                             
	for rows in entries:
		fileLines.append(COLUMN_DELIMITER.join(rows))
		mutalizer.write("".join([rows[1].encode("UTF-8"), "\n"]))
	f.write(ROW_DELIMITER.join(fileLines).encode("UTF-8"))
