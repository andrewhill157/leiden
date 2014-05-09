.. _overview:

Overview
========

About
^^^^^
The Leiden Open Variation Database platform is a popular genetic variant database platform whose installations are home to
many voluntarily curated mutations implicated in a variety of disease areas.

* A list of all current installations: http://www.lovd.nl/2.0/index_list.php
* LOVD platform: http://www.lovd.nl/3.0/home

Unfortunately these variants are in HGVS format (popular in clinical settings) and in coordinates relative to specific
cDNA transcripts, which makes further analysis difficult informatically. Perhaps more concerning is that the standard
for submission of disease causing mutations was has become much stricter in the time since LOVDs inception.
This implies that there are many false positives within this data set. Curation of these databases is completely voluntary,
making many variants completely unreliable or unusable. Despite these challenges, there are likely many true positives
amongst the noise, many of which may not be in other variant databases. Locating and reporting these true positives is
an important goal for the research community.

While LOVD is public access and has provided reST APIs for querying for variants at specific genomic positions or
retrieving some information about LOVD variants in specific genes, none of the available services allow the actual
variant descriptions (or other submitted information) to be downloaded. This package fills that gap and also facilitates
some degree of validation of the data.

The goals of this project are to provide tools for:

* Extracting variants from these databases
* Remapping these variants to VCF format
* Cross-checking of information about these variants to infer concordance of submissions


.. important::
    Actually implicating variants as being pathogenic requires thorough manual curation by examining the full set of information
    (including, but not limited to publication references) for validated variants.  "Validation" as described here
    simply implies correctness and consistency of submitted variants, it does not prove true positive implication in any disease.

General Workflow
^^^^^^^^^^^^^^^^
In general, the workflow is as follows:

1. Extract raw variants from LOVD, saving one tab-delimited file per gene.
2. Annotate variants with VEP (must have VEP on path) and combine with original data in a single VCF file per gene.
3. Validate annotated variants by cross-checking submitted data with annotation and output a single VCF for all variants.

The scripts I have included (see ::ref::driver_scripts_ and ::ref::other_scripts_) make it easy to carry out this workflow.
Custom scripts can also be written for modified workflows using the functionality in the leiden package.

Project Structure
^^^^^^^^^^^^^^^^^

* /bin/ - contains all scripts (see :ref:`driver_scripts` and :ref:`other_scripts` for more info)
* /leiden/ - custom python packages for project (see :ref:`packages` for more info)
* /docs/ - project documentation written in the sphinx framework
* /data/ - sample data formats (see :ref:`data` section for more info)

Other folders are for build and distribution purposes.




