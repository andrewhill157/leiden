.. _data:

Data
====

If you install from source or clone from GitHub, there will be examples of the major data formats used in this project;
variants from the LOVD muscular dystrophy pages.

I have (at the time of this writing) placed a copy of all data extracted from the LOVD muscular dystrophy pages at: http://www.broadinstitute.org/~ahill/lovd_muscle/

Extracted Data
^^^^^^^^^^^^^^
The naming convention is <gene_name>.txt.

Contains original variant data as found on LOVD. Files are saved per gene by extract_data.py.

Raw VCF Files
^^^^^^^^^^^^^
Naming convention is <gene_name>.vcf.

Contains entries from extracted data that successfully remapped to VCF format. Files are save per input file by generate_vcf_files.py.

Annotated VCF Files
^^^^^^^^^^^^^^^^^^^
Naming convention is <gene_name>_ANNOTATED.vcf.

Contains entries from original VCF files with annotations added. Files are saved per input file annotate_vcfs.py.

Final VCF Files
^^^^^^^^^^^^^^^
File is always named validated.vcf and contains variants from all input files.

This is the file format output by validate_annotated_vcf.py. It is a VCF file that contains only validated (in terms of concordance)
variants.

Log Files
^^^^^^^^^
Log files containing information about errors are also saved by some scripts.

remapping_errors.log
++++++++++++++++++++
Contained variants that failed to remap to VCF format and any information about why they failed.

validation_errors.log
+++++++++++++++++++++
Contains variants that caused errors during validation and any information about why they failed.

discordant_annotations.log
++++++++++++++++++++++++++
Contains VCF formatted entries (with annotations) that were found to be discordant during validation.
