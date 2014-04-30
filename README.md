===============================
Leiden Open Variation Database (LOVD) Cleanup
===============================

Tools for extracting, remapping, and validating variants from Leiden Open Variation Database (LOVD) installations.

# Overview of Project Structure

# Overview of General Workflow

# Project Packages

This project contains a number of distinct packages that are responsible for the majority of all functionality in the
project. This does not include any scripts in the root directory.

Note that there are tests for most packages that are not simple wrappers for other scripts or libraries. These tests are
included alongside modules in files called ```test_<module_name>.py```. These tests are not 100% comprehensive, but they
are a good start.

Tests are compatible with the nose unit testing platform. This is an extension of the default unittest platform that makes
tests very easy to develop and run. To run all unit tests for this project run the following from the root directory.
```
nosetests
```

Tests for specific packages/modules can also be run. Please see nose test documentation for more information.

## /lovd

### /lovd/leiden_database.py:

These classes allow a user to extract tables of data (mutations listed for a specific gene in the database) and other
useful information from any Leiden Open Variation Database (LOVD) installation, such as http://www.dmd.nl/nmdb2/. Unfortunately,
it has been necessary to do this by downloading the HTML for relevant pages on the database and parsing out the necessary
data, as they do not provide an easy way to access the data otherwise. Therefore, I have added an external dependency -
beautifulsoup4 - for HTML parsing.

The usage for these classes is as follows:
```
leiden_url = 'http://www.dmd.nl/nmdb2/'  # base URL of LOVD installation
gene_id = 'ACTA1'  # External name for Gene

database = make_leiden_database(leiden_url)  # factory method automatically chooses right database version
database.set_gene_id(gene_id)  # set_gene_id must be called before using the database object

# Get data about the gene
column_labels = leiden_database.get_table_headers()
table_entries = leiden_database.get_table_data()
```
Note that make_leiden_database returns a LeidenDatabase object. There are two subclasses of this type (one for each
version of LOVD). This method ensures that the correct subclass is chosen for the provided URL.

Unit tests for this module are currently quite slow because they actually make requests for HTML data. Ideally, this
would be replaced with Mock Objects where we have canned HTML responses saved on disk.

### /lovd/utilities.py:

These are general utility functions, some of which are used in leiden_database.py.

## /remapping

### macarthur_core/remapping/remapping.py

Genetic variants are often listed in one of two formats, HGVS or VCF. HGVS is compact and has its own (relatively complex)
syntax for describing mutations. However, for large scale analysis projects HGVS is extremely difficult to use effectively
for a number of reasons. This module facilitates conversion back and forth between the two formats.

The class ```VariantRemapper``` in ```/remapping/remapping.py``` wraps a few functions from a third party module
(hgvs - downloadable from PyPI) to make it easier to use within this project. The third party module documentation and
description HGVS vs. VCF notation is described here: https://github.com/counsyl/hgvs.

Unfortunately, the third-party tool depends on two relatively large files that I cannot easily host on github.
These are normally housed in the folder ```/remapping/resources/```. One is a human genome reference sequence (```remapping/resources/hg19.fa```)
and the other is a file containing definitions of transcript sequences that are needed to facilitate conversion between HGVS and VCF
(```/remapping/resources/genes.refSeq```). These two files are hosted at: http://www.broadinstitute.org/~ahill.
These files will need to decompressed using gunzip and placed in ```/remapping/resources/```. The first time these functions
are used, two additional files will be generated (takes some time). Subsequent runs will not require this process to be repeated.

## /input_output

### input_output/file_io.py
This module has functions for reading and writing delimited files to and from 2D lists where first dimension is rows and
the second is columns. It also contains a function for formatting output text in VCF format from list inputs.

### /input_output/web_io.py
This module has functions for reading HTML data from URLs. Essentially, this is just wrapping the library used to
make HTML requests to make it easier to change if needed.

## /broad_cluster
This package contains functions that are useful in the context of running software specifically on the Broad Cluster.

### /broad_cluster/lsf_commands.py
This module contains a wrapper for the lsf command bsub with a number of default parameters. Most (if not all) parameters
can be overridden to customize to user needs.

### /broad_cluster/annotate_vcf.py
This module contains functions that run annotation on VCF files. These include VEP and VCF-Annotate using HGMD/26K/DBSNP.
These are essentially wrappers to call scripts on the broad cluster from the command-line.

## /validation

### /validation/validation.py
This module contains functions that are used to help with validation of the annotated VCF files produced by a driver
script like run_all.py that formats data from LOVD as a VCF and runs specific annotations. These functions rely heavily
on the presence of specific tags in the INFO field of the VCF.

# Driver Scripts

The packages described above provide the main functionality of the project. However, I have also developed a set of
scripts that are specific to my project needs. These simply act as driver scripts to provide a command-line interface
for the protocols I used for this project.

Note that all scripts are implemented using argparse and have built-in help, which accessible via:
```
python <script_name>.py --help
```

## /run_all.py
Run all is the driver script that replicates exactly what I did to produce the set of annotated VCF files that are ready
for validation. Note that because annotation is performed via bsub on the Broad Cluster, the actual validation is not run
by this script. This decision was made because I have very low priority on the cluster, so it was not feasible to block
execution until the annotation step was complete. Future improvements could seek to eliminate this intermediate step.

There are two use-cases for run_all.py:
1. You are starting completely from scratch (no data has been downloaded from LOVD)
```
python run_all.py -u http://www.dmd.nl/nmdb2/ -output_directory my_output_directory
```
This will download data from all genes on the specifed LOVD URL, saving one .txt file (```<gene_name>.txt```) with raw data as
well as two VCF files per gene - one that contains remapped variants in VCF format along with data from LOVD as tags in
the INFO field (```<gene_name>.vcf```) and one that has been annotated with HGMD/26K/DBSNP (```<gene_name>_ANNOTATED.vcf```).

Note that files are not saved for genes with no entries on LOVD. Variants that fail to remap to VCF format are save along
with any information about their failure in ```remapping_errors.log```.

2. You already have the txt files containing raw data from LOVD, but want to re-run the rest of the process. Note that
this was primarily useful during development, but may still have some utility for others.
```
python run_all.py --no_download -output_directory my_output_directory
```

Note that this assumes that the .txt files containing data extracted from LOVD are located in the specified output directory.


# Other Scripts

These scripts are called in sequence by ```run_all.py```. However, they are also callable individually either for debugging
needs or for writing new protocols. Note that the python packages included with this project can also be used to write
entirely scripts if needed. The key point here is that you are not limited to what I have developed here, feel free to
expand and write your own software!

Note that these scripts are listed in the order they are used in my workflow.

## extract_data.py
This is a script that I use in my overall project that makes use of the macarthur_core/lovd and macarthur_core/web_io to extract data from all data from a given LOVD URL. 

The script makes use of argparse to provide a user interface. The string provided for the command-line interface should hopefully provide an explanation of how it is used. It should save a tab-delimited file for each gene's table data, where each output file is named according to the gene name.

From command line generally it is used in the following way:
```
python extract_data.py --all --leiden_url http://www.dmd.nl/nmdb2/ --output_directory results
```
Users can also print a list of all available genes using:
```
python extract_data.py --genes_available --leiden_url http://www.dmd.nl/nmdb2/
```

## generate_vcf_files.py

## annotate_vcfs.py

## validate_annotated_vcfs.py

