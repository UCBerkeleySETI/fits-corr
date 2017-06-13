FitsCorr
========

|License: MIT| |v1.0.0|

FitsCorr is a Python module for aligning and analyzing correlations
between pairs of **ON** and **OFF** fits files.

Table of Contents
-----------------

-  `Installation`_
-  `API`_
-  `Configuration`_
-  `Command-Line`_
-  `Changelog`_

Installation
------------

Install via ``pip`` directly without manually cloning the repository:

.. code:: bash

    $ pip install -e git+https://github.com/UCBerkeleySETI/fits-corr.git#egg=fits-corr

If you'd like to clone the repository, you can run ``pip install`` from within
the cloned directory’s root:

.. code:: bash

    fits-corr $ pip install -e .

If you prefer not to install the package and thus available to your
``$PATH``, the you can build the package as follows:

.. code:: bash

    fits-corr $ python -m fits_corr

API
---

The ``fits_corr`` module exposes a rather simple API. Just point the
correlator to a directory containing the ``.fits`` files or pass in an
**ON** and **OFF** pair instead. You could also specify an output
directory or use the default (``./out/``).

Optionally, you could also enable additional functionalities, which
include enabling coalescing all **JSON** objects to a single file for
easy uploading to a database and multiprocessing for batch processing
the pairs found in the given ``data_dir`` in parallel:

.. code:: python

    from fits_corr import Correlator

    engine = Correlator() \
        .with_data_dir('./your-data-dir/') \
        .with_out_dir('./your-output-dir/') \
        .enable_multiprocessing() \
        .enable_coalescing() \

    engine.run()

Alternatively, you could also supply a **JSON** config file to run the
correlator:

.. code:: python

    from fits_corr import Correlator

    engine = Correlator().load_config('./path-to-config-file')

    engine.run()
    
Configuration
-------------

As FitsCorr supports setup via a **JSON** configuration file, the file
must have the following schema and structure:

.. code:: js

    {
        "coalesce": ( true |false ),
        "align_window": Integer,
        "data_dir": ( "PATH_TO_DATA_DIR" | null ),
        "out_dir": ( "PATH_TO_OUTPUT_DIR" ),
        "histogram_bins": Integer,
        "prmi_radius": Integer,
        "multi_process": ( true | false ),
        "data_pair": ( [ "PATH_TO_ON_PAIR", "PATH_TO_OFF_PAIR" ] | null )
    }

**Note:** Not all of the key-value pairs in the above configuration are
required during setup.

Command-Line
------------

FitsCorr can also be used from the command-line with help and
information available via the ``-h`` flag as follows:

.. code:: bash

    # If installed by pip
    $ fits_corr -h

    # If not installed via pip
    fits-corr $ python fits_corr -h

which yields the following output:

.. code:: text

    usage: fits_corr [-h] (-c DIR | -d DIR | -f ON OFF) [-m] [-g] [-o DIR]
                 [-w WIDTH] [-prmi RADIUS] [-bins BINS] [-s DIR]

    ┌─ Correlation & Similarity Analysis for FITS data files.
    ├─ Version: 0.1.1-alpha
    └─ © Pragaash Ponnusamy 2017

    optional arguments:
      -h, --help            show this help message and exit
      -c DIR, --conf DIR    directory with a _CONFIG.json file to launch with
      -d DIR, --dir DIR     data directory containing .fits files
      -f ON OFF, --file ON OFF
                            path to on and off file pairs

    options:
      -m                    enable multiprocessing
      -g                    enable coalescing
      -o DIR, --out DIR     path to output directory
      -w WIDTH              alignment window width
      -prmi RADIUS          neighborhood window radius for PRMI
      -bins BINS            histogram bins for mutual information
      -s DIR, --save DIR    save directory for config file

Changelog
---------

**Version 1.0.0***

- Minor bug fixes.

**Version 0.1.1-alpha**

- Minor bug fixes.
- Updated usage of list comprehensions.
- Preserved convention of loading configuration file.

**Version 0.1.0-alpha**

-  Initial release.
-  Multiprocessing support for data directory option.
-  Configuration file support.
-  Use of Hermitian matrix for fast SVD.
-  Fast 2d histogram implementation.

.. _Installation: #installation
.. _API: #api
.. _Configuration: #configuration
.. _Command-Line: #command-line
.. _Changelog: #changelog

.. |License: MIT| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
.. |v1.0.0| image:: https://img.shields.io/badge/release-v1.0.0-brightgreen.svg
