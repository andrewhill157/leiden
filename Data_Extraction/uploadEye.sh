# Move results to appropriate folders
mv *.txt LOVD_Eye_Extracted_Data/
mv *.vcf ../Validation/LOVD_Eye_VEP_Validation/

# Transfer VCF files to server
scp ../Validation/LOVD_Eye_VEP_Validation/*.vcf tin:/humgen/atgu1/fs03/ahill/Leiden_Database_Cleanup/eye/
