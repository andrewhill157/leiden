# Move results to appropriate folders
mv *.txt LOVD_Muscular_Dystrophy_Pages_Extracted_Data/
mv *.vcf ../Validation/LOVD_Muscular_Dystrophy_Pages_VEP_Validation/

# Transfer VCF files to server
scp ../Validation/LOVD_Muscular_Dystrophy_Pages_VEP_Validation/*.vcf tin:/humgen/atgu1/fs03/ahill/Leiden_Database_Cleanup/muscle/
