# Move results to appropriate folders
mv *.txt LOVD_Muscular_Dystrophy_Pages_Extracted_Data/
mv *.vcf ../Validation/muscle/

# Transfer VCF files to server
scp ../Validation/muscle/*.vcf tin:/humgen/atgu1/fs03/ahill/Leiden_Database_Cleanup/muscle/
