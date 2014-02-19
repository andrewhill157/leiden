# Move results to appropriate folders
mv *.txt muscle/
mv *.vcf ../validation/muscle/

# Transfer VCF files to server
ssh tin 'rm /humgen/atgu`/fs03/ahill/leiden/muscle/*'
scp ../validation/muscle/*.vcf tin:/humgen/atgu1/fs03/ahill/leiden/muscle/
