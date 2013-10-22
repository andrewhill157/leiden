from extractData import *

# Constant for file type
FILE_EXTENSION = 'txt'
ROW_DELIMITER = '\n'
COLUMN_DELIMITER = ','

# define the gene you want to parse
geneID = 'DYSF'
refSeqID = 'NM_003494.3'

database = LeidenDatabase(geneID, refSeqID)
geneName = database.getGeneName()
fileName = ".".join([geneName, FILE_EXTENSION])

# write table data to file in Unicode encoding (some characters are no ASCII encodable)
with open(fileName, 'w') as f:
		entries = database.getTableData()
		fileLines = []
		headers = database.getTableHeaders()
		fileLines.append(COLUMN_DELIMITER.join(headers))
                                 
		for rows in entries:
                    fileLines.append(COLUMN_DELIMITER.join(rows))
                    
                f.write(ROW_DELIMITER.join(fileLines).encode("UTF-8"))
