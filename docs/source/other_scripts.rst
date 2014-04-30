.. _other_scripts:

Other Scripts
=============

These scripts are called in sequence by ``run_all.py``. However, they are also callable individually either for debugging
needs or for writing new protocols. Note that the python packages included with this project can also be used to write
entirely scripts if needed. The key point here is that you are not limited to what I have developed here, feel free to
expand and write your own software!

Note that these scripts are listed in the order they are used in my workflow.

.. tip::
    Note all scripts are made with argsparse, so contain built-in help. To access help simply execute: ``python <script_name>.py --help``

.. warning::
    These scripts are located alongside the leiden package that contains the packages in this project. In order to use
    packages in scripts elsewhere, the leiden package must be on your PYTHONPATH.

extract_data.py
^^^^^^^^^^^^^^^
extract_data.py allows raw data from the any leiden open variation database installation to be downloaded
and saved to text files (one per gene). There are options to allow either data from all genes or a specific list of genes
to be downloaded as needed. It also allows users to print a list of all available genes at a given URL, which is useful
if you want to check what is available.

.. note::
    Note that both LOVD versions 2 and 3 are supported for this script.

Example Usage
-------------

Download data for all genes from a given url to a specified output directory:

.. code-block:: bash

    python extract_data.py --all --leiden_url http://www.dmd.nl/nmdb2/ --output_directory my_directory

Download a list of specified genes from a given url to a specified output directory:

.. code-block:: bash

    python extract_data.py --leiden_url http://www.dmd.nl/nmdb2/ --output_directory my_directory --gene_list ACTA1 DYSF

Print a list of available genes at a specified URL:

.. code-block:: bash

    python extract_data.py --genes_available --leiden_url http://www.dmd.nl/nmdb2/

generate_vcf_files.py
^^^^^^^^^^^^^^^^^^^^^

generate_vcf_files.py allows the user to generate vcf files for each .txt file output from extract_data.py. Essentially
this file remaps variants from the raw data files to VCF format and then saves VCF formatted files as output.

Note that this script depends on the specific column header names from LOVD 2 installations. Data from these columns is
used both in the remapping process and in creating tags within the INFO column of the VCF. This script could, however,
be adapted to other file formats as needed.

Example Usage
-------------

.. code-block:: bash

    python generate_vcf_files.py -f list_of_file_names.list -o my_output_directory


annotate_vcfs.py
^^^^^^^^^^^^^^^^

.. warning::
    This script depends on resources on the Broad Institute distributed computing cluster, making it somewhat fragile
    and only usable in this context. This decision was made because some of the annotation tools and resources are very large
    and inconvenient to install elsewhere.


annotate_vcfs.py runs annotation on the VCFs produced by generate VCF files. This currently performs the following annotations:
- VEP
- 26K
- HGMD Site, Mutation, Gene
- dbSNP

This script can be adjusted to include additional annotations if needed and can generally be applied to any VCF file.

Example Usage
-------------

.. code-block:: bash

    python annotate_vcfs.py -f input_vcf.vcf -o output_vcf.vcf

.. note::
    This script was made to take single files as input rather than a list because it is useful when submitting jobs
    to lsf on the cluster. Some annotations take quite a while to run.
