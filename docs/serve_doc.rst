Local web interface
===================

``mutagene serve`` starts the local analysis server.
Upload a MAF or VCF file, follow progress, and inspect results in your browser.
The server uses the same analysis code as the command-line subcommands.

Install and start
-----------------

The web interface requires additional dependencies:

.. code-block:: bash

   pip install mutagene[web]

Without these dependencies, ``mutagene serve`` reports the missing packages and exits.

.. warning::

   The server has no authentication.
   Keep the default ``127.0.0.1`` address unless you control access from other machines.

Start the server:

.. code-block:: bash

   mutagene serve

Use another port without opening a browser:

.. code-block:: bash

   mutagene serve --port 8080 --no-browser

Select the reference genome
---------------------------

Choose the assembly that matches your mutation coordinates.
The server stores assemblies in ``~/.mutagene/genomes/`` as ``<name>.2bit``.
Download an assembly before analysis:

.. code-block:: bash

   mutagene fetch genome hg19

The interface lists local assemblies and can download a missing assembly on request.
The hg19 and hg38 files each require roughly 800 MB.
The analysis compares input reference alleles with the assembly.
A high mismatch rate triggers a warning on the results page.

Server options
--------------

``--host HOST``
    Set the listening address. The default is ``127.0.0.1``.
    ``0.0.0.0`` allows connections from other machines.

``--port PORT``
    Set the port. The default is ``5000``.

``--no-browser``
    Start without opening a browser window.

``--debug``
    Use only during development.
    This option exposes Flask's interactive debugger to anyone who can reach the server.

``--params-out FILE``
    Save run parameters as JSON.

``--params-in FILE``
    Read parameters from a file that ``--params-out`` creates.

Find stored results
-------------------

==========================  =================================================
Path                        Contents
==========================  =================================================
``~/.mutagene/results.db``  Analyses, their status and their results (SQLite)
``~/.mutagene/results/``    Output files, one directory per analysis
``~/.mutagene/uploads/``    Uploaded input files
``~/.mutagene/genomes/``    Downloaded 2bit assemblies
==========================  =================================================

The history page retains analyses across server restarts.
You can inspect completed results without uploading the input again.

When you repeat an analysis, the server writes new output to a staging directory.
It replaces the stored output only after the new run succeeds.
A failed run preserves the previous results.
After an interruption during replacement, the server reconciles files with the database at startup.
The committed database status determines which output it retains.
