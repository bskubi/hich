Hich CLI Utilities
==================

``hich compartments [options] REFERENCE MATRIX RESOLUTION``
.................................................................

Calls compartment scores on Open2C-format Hi-C data.

Options:

+ ``--chroms``: Chrom/contig names to use. If given, contigs not named here will not be used.
+ ``--exclude-chroms``: Chrom/contig names to exclude. If given, contigs named here will not be used, even if specified in ``--chroms``.
+ ``--keep-chroms-when``: Python code that references a string ``chrom`` and evaluates to True if that chromosome name should be used.
+ ``--n_eigs``: Number of eigenvectors to call. These will be sign-flipped so that positive scores correspond to regions of higher %GC content to create compartment scores. Default: 1

Arguments:

+ ``REFERENCE``: Genome fasta reference used to determine %GC content.
+ ``MATRIX``: A .cool or .mcool file used as input to `cooltools eigs-cis <https://cooltools.readthedocs.io/en/latest/cli.html#cooltools-eigs-cis>`_
+ ``RESOLUTION``: Resolution at which compartment scores are called.

Example:

.. code-block:: bash
    
    hich compartments --keep-chroms-when "chrom.startswith('chr')" hg38_noalts.fa.gz k562.mcool 10000

``hich digest [options] REFERENCE DIGEST``
.........................

Creates restriction enzyme fragment index in BED format. 

Options:

+ ``--output``: Output file. Compression autodetected by file extension. If None, prints to stdout.
+ ``--startshift``: Fixed distance to shift start of each fragment. Default: 0
+ ``--endshift``: Fixed distance to shift end of each fragment. Default: 0
+ ``--cutshift``: Fixed distance to shift end of each fragment. Default: 1

Arguments:

+ ``REFERENCE``: Reference genome to use as basis for digest
+ ``DIGEST``: Space-delimited list of enzyme names and/or kit name

Examples:

.. code-block:: bash
    
    hich digest --output hg38_arima.bed hg38_noalts.fa.gz Arima
    hich digest --output hg38_hic3.bed hg38_noalts.fa.gz DpnII DdeI

Supported enzymes:
+ All of ~800 REBASE enzymes used by `biopython's restriction module <https://github.com/biopython/biopython/blob/master/Doc/cookbook/Restriction/Restriction.md>`_, as well as ``Arima``, ``Phase Plant``, ``Phase Animal``, ``Phase Microbiome``, ``Phase Human``, ``Phase Fungal``.

``hich downsample [options] ``
.........................

Options:

+ ``--conjuncts``:
+ ``--cis-strata``:
+ ``--orig-stats``:
+ ``--target-stats``:
+ ``--to-size``:

Arguments:

+ ``INPUT_PAIRS_PATH``
+ ``OUTPUT_PAIRS_PATH``


``hich fragtag``
.........................

``hich gather``
.........................

``hich hicrep``
.........................

``hich organize``
.........................

``hich reshape``
.........................

``hich stats_aggregate``
.........................

``hich stats``
.........................

