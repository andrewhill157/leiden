.. _packages:

Packages
========
This project contains a number of distinct packages that are responsible for the majority of all functionality in the
project. This does not include any scripts in the root directory.

Note that there are tests for most packages that are not simple wrappers for other scripts or libraries. These tests are
included alongside modules in files called ``test_<module_name>.py``. These tests are not 100% comprehensive, but they
are a good start.

Tests are compatible with the nose unit testing platform. This is an extension of the default unittest platform that makes
tests very easy to develop and run. To run all unit tests for this project run the following from the root directory.

.. code-block:: python

    nosetests


Tests for specific packages/modules can also be run. Please see nose documentation for more information.

lovd
^^^^

leiden_database
---------------
These classes allow a user to extract tables of data (mutations listed for a specific gene in the database) and other
useful information from any Leiden Open Variation Database (LOVD) installation, such as http://www.dmd.nl/nmdb2/. Unfortunately,
it has been necessary to do this by downloading the HTML for relevant pages on the database and parsing out the necessary
data, as they do not provide an easy way to access the data otherwise. Note that I have chosen to use beautifulsoup4 internally
for HTML parsing.

The usage for these classes is as follows:

.. code-block:: python

    leiden_url = 'http://www.dmd.nl/nmdb2/'  # base URL of LOVD installation
    gene_id = 'ACTA1'  # External name for Gene

    database = make_leiden_database(leiden_url)  # factory method automatically chooses right database version
    database.set_gene_id(gene_id)  # set_gene_id must be called before using the database object

    # Get data about the gene
    column_labels = leiden_database.get_table_headers()
    table_entries = leiden_database.get_table_data()

Note that make_leiden_database returns a LeidenDatabase object. There are two subclasses of this type (one for each
version of LOVD). This method ensures that the correct subclass is chosen for the provided URL.

Unit tests for this module are currently quite slow because they actually make requests for HTML data. Ideally, this
would be replaced with Mock Objects where we have canned HTML responses saved on disk.

Member Descriptions
+++++++++++++++++++

.. automodule:: leiden.lovd.leiden_database
   :members:


utilities
---------
These are general utility functions, some of which are used in leiden_database.py as general purpose functions.

Member Descriptions
+++++++++++++++++++

.. automodule:: leiden.lovd.utilities
   :members:

remapping
^^^^^^^^^
Genetic variants are often listed in one of two formats, HGVS or VCF. HGVS is compact and has its own (relatively complex)
syntax for describing mutations. However, for large scale analysis projects HGVS is extremely difficult to use effectively
for a number of reasons. This module facilitates conversion back and forth between the two formats.

The class ``VariantRemapper`` in ``/remapping/remapping.py`` wraps a few functions from a third party module
(hgvs - downloadable from PyPI) to make it easier to use within this project. The third party module documentation and
description HGVS vs. VCF notation is described here: https://github.com/counsyl/hgvs.

.. important::
    Unfortunately, the third-party tool used for remapping depends on a relatively large file that I cannot easily host on Github.
    This is normally housed in the folder ``/remapping/resources/``. It is a human genome reference sequence (``remapping/resources/hg19.fa``)
    I have temporarily hosted a copy at at: http://www.broadinstitute.org/~ahill. This file will need to decompressed using gunzip
    and placed in ``/remapping/resources/``. The first time this package is used, two additional files will be generated
    (takes some time). Subsequent runs will not require this process to be repeated.

remapping
---------

Member Descriptions
+++++++++++++++++++

.. automodule:: leiden.remapping.remapping
   :members:

input_output
^^^^^^^^^^^^

file_io
-------
This module has functions for reading and writing delimited files to and from 2D lists where first dimension is rows and
the second is columns. It also contains a function for formatting output text in VCF format from list inputs.

Member Descriptions
+++++++++++++++++++

.. automodule:: leiden.input_output.file_io
   :members:

web_io
------
This module has functions for reading HTML data from URLs. Essentially, this is just wrapping the library used to
make HTML requests to make it easier to change if needed.

Member Descriptions
+++++++++++++++++++

.. automodule:: leiden.input_output.web_io
   :members:


broad_cluster
^^^^^^^^^^^^^
This package contains functions that are useful in the context of running software specifically on the Broad Institute
distributed computing cluster.

.. warning::
    These packages depend on resources on the Broad Institute distributed computing cluster, making them somewhat fragile
    and only usable in this context. This decision was made because some of the annotation tools and resources are very large
    and inconvenient to install elsewhere.

lsf_commands
------------
This module contains a wrapper for the lsf command bsub with a number of default parameters. Most (if not all) parameters
can be overridden to customize to user needs.

Member Descriptions
+++++++++++++++++++

.. automodule:: leiden.broad_cluster.lsf_commands
   :members:

annotate_vcf
------------
This module contains functions that run annotation on VCF files. These include VEP and VCF-Annotate using HGMD/26K/DBSNP.
These are essentially wrappers to call scripts on the broad cluster from the command-line.

Member Descriptions
+++++++++++++++++++

.. automodule:: leiden.broad_cluster.annotate_vcf
   :members:

validation
^^^^^^^^^^

validation
----------
This module contains functions that are used to help with validation of the annotated VCF files produced by a driver
script like run_all.py that formats data from LOVD as a VCF and runs specific annotations. These functions rely heavily
on the presence of specific tags in the INFO field of the VCF.

.. automodule:: leiden.validation.validation
    :members: