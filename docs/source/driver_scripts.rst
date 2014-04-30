.. _driver_scripts:

Driver Scripts
==============

The packages described above provide the main functionality of the project. However, I have also developed a set of
scripts that are specific to my project needs. These simply act as driver scripts to provide a command-line interface
for the protocols I used for this project.

Note that all scripts are implemented using argparse and have built-in help, which accessible via:

.. code-block:: bash

    python <script_name>.py --help

run_all.py
^^^^^^^^^^
run_all.py is the driver script that replicates exactly what I used to produce the set of annotated VCF files that are ready
for validation.

.. warning::
    This script depends on resources on the Broad Institute distributed computing cluster, making it somewhat fragile
    and only usable in this context. This decision was made because some of the annotation tools and resources are very large
    and inconvenient to install elsewhere.

.. important::
    Note that because annotation is performed via bsub on the Broad Institute distributed computing cluster, the actual validation is not run
    by this script. This decision was made because I have very low priority on the cluster, so it was not feasible to block
    execution until the annotation step was complete. Future improvements could seek to eliminate this intermediate step.

Example Usage
-------------

There are two use-cases for run_all.py:

1. You are starting completely from scratch (no data has been downloaded from LOVD)

.. code-block:: bash

    python run_all.py -u http://www.dmd.nl/nmdb2/ -output_directory my_output_directory

This will download data from all genes on the specifed LOVD URL, saving one .txt file (``<gene_name>.txt``) with raw data as
well as two VCF files per gene - one that contains remapped variants in VCF format along with data from LOVD as tags in
the INFO field (``<gene_name>.vcf``) and one that has been annotated with HGMD/26K/DBSNP (``<gene_name>_ANNOTATED.vcf``).

Note that files are not saved for genes with no entries on LOVD. Variants that fail to remap to VCF format are save along
with any information about their failure in ``remapping_errors.log``.

2. You already have the txt files containing raw data from LOVD, but want to re-run the rest of the process. Note that
this was primarily useful during development, but may still have some utility for others.

.. code-block:: bash

    python run_all.py --no_download -output_directory my_output_directory


Note that this assumes that the .txt files containing data extracted from LOVD are located in the specified output directory.

validate_annotated_vcfs.py
^^^^^^^^^^^^^^^^^^^^^^^^^^
validate_annotated_vcfs.py takes a list of annotated VCFs (such as those output by run_all.py or annotate_vcfs.py (see :ref:`other_scripts`)
and produces a single VCF file containing all variants that were confirmed to be concordant. This script tries to cross
reference information from the annotation with information provided in the entries from LOVD.

All variants are saved to a single output file as specified by the user.

* Enough information must be provided to perform validation
* Mutations that cause amino acid change: REF and ALT much match annotation
* Synonymous: annotation must also predict synonymous mutations.
* Other:
    - Splice: Mutation must match annotation. Must also be in conserved splice region at beginning or end of exon and match expectation based on conserved patterns.
    - Frameshift: ensembl API is used to confirm original amino acid, alternate amino acid, and location of any stop codons.
    - Codon Loss: ensembl API is used to confirm original amino acid and alternate amino acid.

Example Usage
-------------

.. code-block:: bash

    python validate_annotated_vcfs.py -f input_file_list.list -o output_file.vcf