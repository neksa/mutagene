=============================================
Serve: the local web interface
=============================================

-----------------
1. Description
-----------------

``mutagene serve`` starts a small web application on your own machine for
running analyses interactively: upload a MAF or VCF file, watch the analysis
progress, and browse the resulting profiles, signature decompositions,
clustering and driver rankings in the browser.

It runs the same analysis code as the command-line subcommands rather than a
separate implementation, so results agree with what ``profile``, ``signature``
and ``rank`` produce for the same input.

The server binds to ``127.0.0.1`` by default, so it is reachable only from your
own machine. It has no authentication; do not expose it to a network you do not
control.

-----------------
2. Installation
-----------------

The web interface needs extra dependencies that the base package does not
install::

    pip install mutagene[web]

Without them, ``mutagene serve`` reports what is missing and exits rather than
failing with an import error.

-----------------
3. Arguments
-----------------

**Command:** ``mutagene serve [arguments]``

--host HOST
    Address to bind to. Default ``127.0.0.1``. Use ``0.0.0.0`` only if you
    intend the server to be reachable from other machines.

--port PORT
    Port to bind to. Default ``5000``.

--no-browser
    Do not open a browser window on startup.

--debug
    Enable Flask debug mode. Development only: it exposes an interactive
    debugger to anyone who can reach the server.

--params-out FILE
    Write the parameters of this run to FILE as JSON.

--params-in FILE
    Take parameters from a FILE written by ``--params-out``.

Examples::

    mutagene serve
    mutagene serve --port 8080 --no-browser

-----------------
4. Genomes
-----------------

An analysis needs the 2bit assembly its coordinates refer to. Assemblies live in
``~/.mutagene/genomes/`` as ``<name>.2bit`` and can be fetched in advance::

    mutagene fetch genome hg19

The web interface offers the assemblies it finds there and can download a
missing one on request. The files are large -- roughly 800 MB each for hg19 and
hg38 -- so the first download takes a while.

Choosing the wrong assembly is the most common cause of implausible results. The
analysis compares each reference allele against the assembly and shows a warning
in the results page when a large share of them disagree, naming the assembly to
try instead.

-------------------------------
5. Where results are stored
-------------------------------

==========================  =================================================
Path                        Contents
==========================  =================================================
``~/.mutagene/results.db``  Analyses, their status and their results (SQLite)
``~/.mutagene/results/``    Output files, one directory per analysis
``~/.mutagene/uploads/``    Uploaded input files
``~/.mutagene/genomes/``    Downloaded 2bit assemblies
==========================  =================================================

Analyses persist across restarts and remain listed in the history page, so a
completed analysis can be revisited without re-uploading its input.

Re-running an analysis replaces its stored output only once the new run has
succeeded: the results are written aside and swapped into place together with
their database rows. A failed re-run therefore leaves the previous results
intact and displayable, and a server that is interrupted mid-swap reconciles the
files with the database on next startup.
