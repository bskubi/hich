Configuration
=============

We will distinguish between **Nextflow config options (NF-config)**, **Hich config options (Hich-config)**, and **sample attributes**. NF-config controls the Nextflow workflow management system (the platform on which Hich runs) and are listed `here <https://www.nextflow.io/docs/latest/reference/config.html>`_. Hich-config controls the way Hich processes individual datasets and builds **samples**, which are bundles of sample attributes including paths to data files containing data used for processing (i.e. fastq files, aligner index files, contact matrices, etc.) as well as directives for how they should be processed. (i.e. minimum MAPQ filters or how loops should be called).

To set up a run, you will typically use three configuration files. As sample attributes can be specified at the command line, none are required. A maximum of one configuration file of each type can be used.

1. To specify sample attributes to individual samples, the easiest way is with a tab-separated value (TSV) **sample file** , specified with the ``--sample-file`` option. Examples in the ``vignettes`` directory.
2. To specify Hich-config offering more powerful control over how Hich assigns attributes to samples, such as global default values or complex parameter exploration plans, you can use a `YAML <https://en.wikipedia.org/wiki/YAML>`_ **params file**, specified with the ``-params-file`` option. Examples in the ``params`` directory.
3. To define NF-config, such as computational resource management profiles and control how outputs are published, as well as a few rarely modified Hich-config options, a `Nextflow configuration syntax <https://www.nextflow.io/docs/latest/config.html>`_ **Nextflow config file**, typically named ``nextflow.config`` or specified at the command line with the ``-c`` option. Example in the ``hich`` directory. 

Hich also allows many configuration options to be specified at the command line, using options like ``--defaults.aligner bwa-mem2``.

The `Nextflow training page on configuration options <https://training.nextflow.io/2.0/nf_customize/04_config/>`_ page describes how conflicts are resolved.

.. note::
    NF-config uses a ``-`` prefix, as for options like ``-params-file`` and ``-c``. Hich-config uses a double ``--`` prefix, as for options like ``--sample-file``.

.. _sample-file:
Sample file
.........................

The sample file is typically a tab-separated value (TSV) file that describes sample-specific attributes.

``--sampleFile``
    The path to the sample file.

Formatting requirements:

+ Delimiter matches ``--sampleFileSep`` parameter
+ Has headers for every column with content. No content-free headerless columns between columns with content.
+ Each row corresponds to one sample.
+ Each column corresponds to one sample attribute, with the header name matching the name expected by Hich for that sample attribute.

.. note::
    You can choose a different delimiter by setting the ``--sampleFileSep`` parameter, like ``--sampleFileSep ","`` to use a comma-separated value (CSV) format for the sample file. This is specified in ``nextflow.config`` by default, and can be modified there or overridden at the command line.

Params file
..................

The params file is in `YAML syntax <https://en.wikipedia.org/wiki/YAML>`_ and offers powerful ways to specify sample attributes and parameterize the workflow:

+ Attribute values applied to all samples, or a subset of samples, by default, by specifying values of sample attributes in its ``defaults`` block
+ Complex parameterizations for individual workflow steps
+ Relationships between samples to control downsampling, merging, and deduplication

Examples of premade params files can be found in the ``hich/params`` directory.

Nextflow config file
.........................

The Nextflow config file is in `Nextflow configuration syntax <https://www.nextflow.io/docs/latest/config.html>`_ and is typically named ``nextflow.config``. It offers exactly the same capabilities as the params file. The reason for using both is for convenience: infrequently changed config options that are shipped with Hich are in ``nextflow.config``, while config options more likely to be adjusted to suit individual runs are placed in params files. The ``nextflow.config`` file included with Hich includes a number of resource management profiles in the ``profiles`` section, as well as more general parameters in the ``params`` and ``params.general`` sections, including which containers are used, where outputs are published, when to generate multiQC reports, how many reads to keep for humid runs, and what separator delimits columns in the sample file.

Set samples at the command line
.........................................

The following allow including fastq samples, and even sample-specific attributes encoded in the fastq filename, directly from the command line. Samples declared this way will be used in addition to the :ref:`sample file <sample-file>`, if specified.

``--fastqInterleaved``
    Interleaved fastq files (i.e. r1 followed by r2 in the same file). Example: ``--fastqInterleaved fq/*.fq.gz``

``--fastqPairs``
    Paired fastq files (i.e. an r1 and an r2 file). Filenames are parsed using Nextflow's `fromFilePairs <https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs>`_ syntax. Example: ``--fastqPairs fq/*.r{1,2}.fq.gz``

``--samplesFromSRA``
    URL to paired fastq files hosted on `SRA <https://www.ncbi.nlm.nih.gov/sra>`_.

``--paramsFromPath``
    Uses a similar syntax to Python's `parse <https://pypi.org/project/parse/>`_ library to extract sample attribute values from filenames using a parsing pattern. Note that if using ``--fastqPairs``, only the first file will be parsed (i.e. r1 if that was specified first between brackets, as in ``*.r{1,2}.fq.gz``). Example: ``--fastqPairs fq/*.r{1,2}.fq.gz --paramsFromPath {condition}_{biorep}_{techrep}.r1.fq.gz`` 

``--samples``
    Interprets filename to read in fastq (``".fastq", ".fq"``), sam/bam (``".sam", ".bam"``), pairs (``".pairs"``), mcool (``".mcool"``), or hic (``".hic"``) files. The extension can be included in the filename (``"*.fq.gz"``) or be at the end of the filename (``*.fq``). Example: ``--samples data/*.pairs.gz``

Set Hich-config at the command line
...........................................

Hich-config, such as that specified in the :ref:`params file <params-file>` or :ref:`nextflow.config <nextflow-config-file>`, can also be specified directly via the command line. This will override that option's specification in the params file or nextflow.config file. Example: ``nextflow run contents/main.nf --defaults.minMapq 10 --general.publish.mode copy``

.. include:: hich_config.rst

.. include:: sample_attributes.rst