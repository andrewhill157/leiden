# Move results to appropriate folders
mv *.txt eye/
mv *.vcf ../validation/eye/

# Transfer VCF files to server
scp ../validation/eye/*.vcf tin:/humgen/atgu1/fs03/ahill/leiden/eye/
