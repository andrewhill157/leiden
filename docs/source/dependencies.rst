.. _dependencies:

Installation
============

Installation for General Use (No Development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While this project is registered in PyPI and ideally could be installed using ``pip install leiden``, two dependencies
(pygr and hgvs) do not always install correctly with pip.

In the next sections I will step through the installation process for the two finicky dependencies.

.. important::
    Install the dependencies below first and then install leiden.

Once these two initial dependencies have been installed, you can do any of the following:

Clone:

.. code-block:: bash

    git clone git@github.com:andrewhill157/leiden.git

Clone and Install from Source:

.. code-block:: bash

    git clone git@github.com:andrewhill157/leiden.git
    cd leiden
    python setup.py install leiden

Install using pip:

.. code-block:: bash

    pip install leiden


All project dependencies are listed in requirements.txt. If you find that you are missing any requirements after the
installation you can install them using:

.. code-block:: bash

    pip install -r requirements.txt

pygr
----

pygr requires that you have bsddb3 and Berkeley DB installed. This can be installed using homebrew:

.. code-block:: bash

    brew install berkeley-db --without-java

    # Note prefix to pip not required on Broad cluster
    sudo BERKELEYDB_DIR=/usr/local/Cellar/berkeley-db/5.3.15/ pip install bsddb3


.. tip::
    bsddb3 and Berkeley DB should already be installed on the Broad cluster!

Once bbd3 and Berkeley DB are installed, you should be able to install pygr with ``pip install pygr``. If there are still errors,
the following may correct the problem.

.. code-block:: bash

   export CFLAGS=-Qunused-arguments
   export CPPFLAGS=-Qunused-arguments
   pip install pygr

I have not encountered any additional problems installing pygr.

hgvs
----

I have had no luck installing hgvs with pip. I am not sure exactly what is causing the problem, so I advise cloning
the Github repository and installing from source:

.. code-block:: bash

    git clone git@github.com:counsyl/hgvs.git
    python hgvs setup.py install
    rm -rf hgvs

I have not encountered any additional problems installing hgvs.

.. important::
    Unfortunately, this tool depends on a relatively large file that I cannot easily host on Github.
    This is normally housed in the folder ``/leiden/remapping/resources/``. It is a human genome reference sequence (``hg19.fa``)
    I have temporarily hosted a copy at at: http://www.broadinstitute.org/~ahill. This file will need to decompressed using gunzip
    and placed in ``/leiden/remapping/resources/``. The first time this package is used, two additional files will be generated
    (takes some time). Subsequent runs will not require this process to be repeated.

Installation for Development
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you would like to extend or modify the existing code-base or scripts while still having the package installed,
you can install in editable or development mode. This differs slightly from the default installation mode.

.. important::
    Just as with installation for general use, pygr and hgvs must be installed prior to running the commands below. See
    Installation for General Use section for more information.

The easiest way to do this is to install from cloned source.

.. code-block:: bash

    git clone git@github.com:andrewhill157/leiden.git
    cd leiden
    python setup.py develop

Either of these methods will make the leiden packages accessible by python, but allow you to edit and call the modified
source without re-installing the package. Note that the dependencies must still be installed according to instructions in
the Installation for General Use section.

If you want to contribute to the project via Github, please see the :ref:`contributing` page.