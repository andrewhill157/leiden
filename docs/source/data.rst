.. _data:

Data
====

If you install from source or clone from GitHub, there will be examples of the major data formats output from the scripts
in this package.

Extracted Data
^^^^^^^^^^^^^^
The naming convention is <gene_name>.txt.

Contains original variant data as found on LOVD. Files are saved per gene by extract_data.py.

Annotated VCF Files
^^^^^^^^^^^^^^^^^^^
Naming convention is <gene_name>_ANNOTATED.vcf.

Contains entries from original VCF files with annotations added. Files are saved per input file annotate_vcfs.py.

Final VCF Files
^^^^^^^^^^^^^^^
This is the file format output by validate_annotated_vcfs.py. It is a VCF file that contains only concordant variants.

Log Files
^^^^^^^^^
Log files containing information about errors are also saved by some scripts.

remapping_errors.log
++++++++++++++++++++
Contained variants that failed to remap to VCF format and any information about why they failed.

discordant_annotations.vcf
++++++++++++++++++++++++++
Contains VCF formatted entries (with annotations) that were found to be discordant during validation.
