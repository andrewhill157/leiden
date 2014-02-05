# Move results to appropriate folders
mv *.txt LOVD_Eye_Extracted_Data/
mv *.vcf ../Validation/eye/

# Transfer VCF files to server
scp ../Validation/eye/*.vcf tin:/humgen/atgu1/fs03/ahill/Leiden_Database_Cleanup/eye/
