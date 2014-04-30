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
cDNA transcripts. Perhaps more concerning is that the standard for submission of disease causing mutations was has become much
stricter in the time since LOVDs inception. This implies that there are many false positives within this data set. Voluntary
curation of these databases is completely voluntary, making many variants completely unreliable or unusable. Despite these
challenges, there are likely many true positives amongst the noise, many of which may not be in HGMD. Locating and reporting
these true positives is an important goal for the research community.

A technical barrier to this effort is that bioinformatics researchers need these variants to be in VCF format because these
coordinates are in genomic space and compatible with popular tools such as VEP. While LOVD Version 3 is currently seeking
to address this problem, many databases have not been migrated. Furthermore, LOVD Version 3 lists variants in HGVS genomic coordinates, not VCF format.
Lastly, much of the variants (or important fields from their entries) are not easily accessible publicly.

The goals of this project are to provide tools for:

* Extracting variants from these databases
* Remapping these variants to VCF format
* Cross-checking of information about these variants to infer concordance of entries
* Augmenting existing data with allele frequency data to allow filtering of common variants

This is a good first step towards finding true variants. Actually implicating variants as being pathogenic will require
manual curation by examining any publication references included with the extracted data. This is a future goal and not
within the scope of the tools provided here.

General Workflow
^^^^^^^^^^^^^^^^
In general, the workflow is as follows:

1. Extract raw variants from LOVD
2. Remap raw variants to VCF format
3. Annotate VCF files
4. Validate annotated variants by cross-checking submitted data with annotation

The scripts I have included make it easy to carry out this workflow. Custom scripts can also be written for modified
workflows using the functionality in the leiden package.

Project Structure
^^^^^^^^^^^^^^^^^

* /bin/ - contains all scripts (see :ref:`driver_scripts` and :ref:`other_scripts` for more info)
* /leiden/ - custom python packages for project (see :ref:`packages` for more info)
* /docs/ - project documentation written in the sphinx framework
* /data/ - sample data formats (see :ref:`data` section for more info)

Other folders are for build and distribution purposes.




