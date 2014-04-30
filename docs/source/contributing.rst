.. _contributing:

Contributing to Leiden
======================

The project is hosted on `Github <https://github.com/andrewhill157/leiden>`_ and can be cloned with the following command:

.. code-block:: bash

    git clone git@github.com:andrewhill157/leiden.git

.. note::
    Even if you do not install the package, you must install in development mode or add the leiden package to your PYTHONPATH
    to use leiden (see :ref:`dependencies` page).

If you fork my project on Github, I will gladly review any pull requests. If you want your changes to be incorporated into
the main project, this is the best approach. Please include a detailed description with any requests.

The lovd and remapping packages and their respective scripts are usable and extensible by anyone!

Unfortunately, the annotation and validation packages and scripts that use them are tightly dependent on resources only
available on the Broad Institute distributed computing cluster or their output format. This decision was made because some of the annotation tools
and resources are very large and inconvenient to install elsewhere. Therefore, only members of the Broad Institute will able to make
use of this particular code. Feedback on how to eliminate this coupling is always appreciated.