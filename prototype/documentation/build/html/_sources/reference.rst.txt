Hich Reference
==============

Sample file
-----------

Normally, the `sample file`_ is called "samples.tsv" (tab-delimited). Basic `sample attributes`_ are usually specified here. Default sample attributes are customizeable and can be specified on the basis of individual sample ids in `nextflow.config`_, which is useful for specifying defaults for `biorep`_ and `condition`_ samples produced via merge.

.. csv-table:: **Example 1.** Because `reference`_, `chromsizes`_, `index_dir`_, `index_prefix`_, and `fragfile`_ files are unspecified and the `assembly`_ values are supported, Hich will download the reference and produce these needed files automatically.
    :file: tables/reference_samplefile1.tsv
    :delim: tab
    :header-rows: 1

.. csv-table:: **Example 2.** Here, needed reference files are given (possibly from a permanent lab repository), so they will be used rather than produced by Hich. Because there's just one sample, there is no need to specify a `biorep`_ or `techrep`_ parameter.
    :file: tables/reference_samplefile2.tsv
    :delim: tab
    :header-rows: 1

.. csv-table:: **Example 3.** Here, files with a mixture of `datatype`_ values are ingested by Hich.
    :file: tables/reference_samplefile3.tsv
    :delim: tab
    :header-rows: 1

.. csv-table:: **Example 4.** An experiment using a variety of `enzymes`_ for reference digestion and fragment tagging, as well as one sample not tagged or filtered (MNase).
    :file: tables/reference_enzymes.tsv
    :delim: tab
    :header-rows: 1

nextflow.config
---------------

The nextflow.config file is one way to `configure <https://www.nextflow.io/docs/latest/config.html>`_ Nextflow, including by setting Hich-specific sample attributes. All sample attributes are described in this section.

Scopes
......

Hich uses specialized config scopes, specified with a name followed by brackets, to group related sample attributes and general Hich workflow parameters. Here is an example with a subset of the real Hich default ``nextflow.config`` file and an extra scope used to specify parameters for a merge.

.. code-block:: c

    params {
        general {
            // The general scope holds params
            // relevant to general Hich workflow
            // control, not sample attributes.

            sampleFile {
                // Path to sample file
                // and column separator.
                filename = "samples.tsv"
                sep = "\t"
            }
        }

        defaults {
            // The default scope gives default
            // sample attributes applied to all
            // samples if an explicit value is not
            // given in samples.csv or in a scope
            // specific to the sample's id.

            // Default techrep and biorep labels
            // applied to any samples where they
            // are not specified in samples.csv
            techrep = 1
            biorep = 1

            // Minimum mapq threshold to keep reads
            min_mapq = 30
        }

        ko {
            // Apply these parameters to samples
            // with the id "KO" and "NT"
            ids = ["KO", "NT"]
            hicrep {
                exclude = ["chrM"]
            }
        }
    }

general
,,,,,,,

last_step
^^^^^^^^^
Specifies the last processing step that should be executed when ``nextflow run hich.nf`` is invoked (as a stub, humid or full run). QC for that step will also be completed. Useful for test runs, debugging, and making processing decisions based on QC results. Commented out by default.

.. code-block:: c

    params {
        general {
            //last_step = "align"

sampleFile
^^^^^^^^^^
The filename and column separator for the `sample file`_. The ``filename`` param can contain a path relative to the Nextflow `projectDir <https://www.nextflow.io/docs/latest/script.html>`_.

.. code-block:: c

    params {
        general {
            sampleFile {
                filename = "samples.tsv"
                sep = "\t"
            }

publish
^^^^^^^
Specifies the Nextflow `publishDir <https://www.nextflow.io/docs/latest/process.html#publishdir>`_ mode and output directory for the results of various Hich processes.

.. code-block:: c

    params {
        general {
            publish {
                // Nextflow publishDir param for all processes
                // https://www.nextflow.io/docs/latest/process.html#publishdir
                mode = "copy"

                // Where to publish results of Hich processes
                genome = "resources/.hich"
                chromsizes = "resources/.hich"
                bwa_mem2_index = "resources/.hich/bwa-mem2/index"
                bwa_mem_index = "resources/.hich/bwa-mem/index"
                digest = "resources/.hich"

                bam = "results/bam"
                parse = "results/pairs/parse"
                dedup = "results/pairs/dedup"
                mcool = "results/matrix/mcool"
                hic = "results/matrix/hic"
                pair_stats = "results/pair_stats"
                qc = "results/qc"
            }

qc_after
^^^^^^^^

.. code-block:: c

    params {
        general {
            // After these steps, generate read-level pairs 
            // stats files and generate a combined MultiQC report
            // for all samples at each processing stage
            qc_after = ["Parse",
                        "IngestPairs",
                        "OptionalFragtag",
                        "TechrepsToBioreps",
                        "Deduplicate",
                        "BiorepsToConditions",
                        "Select"]

humid
^^^^^

.. code-block:: c

    params {
        general {
            // Number of reads to downsample to
            // when doing a humid run
            humid {
                n_reads = 100000
            }

defaults
,,,,,,,,

All `sample attributes`_ specified under this scope will be applied to any samples for which a value is not given in the `sample file`_ or one of the `custom scopes`_.

custom scopes
,,,,,,,,,,,,,

Custom scopes work just like the `defaults`_ scope, except that they have a special ``ids`` list specifying the set of ids to which they should be applied. Custom scopes override the values in the `sample file`_.

Sample attributes
-----------------

In Hich, a **sample** is a single unit of data, such as a **technical replicate**, **biological replicate**, or experimental **condition**. Each sample has a number of sample attributes. These can be specified via columns in the `sample file`_, or to a subset of sample ids via the `nextflow.config`_ file (or anywhere your Nextflow is `configured <https://www.nextflow.io/docs/latest/config.html>`_, including directly at the command line).


Basic
.....
condition
,,,,,,,,,,,,
| Required (no default)
| 
| A label for the condition. Biological replicates with the same `condition`_ label will be merged into a condition sample.

biorep
,,,,,,,,,
| Required (default = ``1``)
| 
| A label for the biological replicate. Technical replicates with the same `condition`_ and `biorep`_ labels will be merged into a biorep sample. Note that Hich does not increment the default value, so it is essential to explicitly specify a biorep label if a value other than 1 is desired.

techrep
,,,,,,,,,,
| Required (default = ``1``)
| 
| A label for the technical replicate. Note that Hich does not increment the default value, so it is essential to explicitly specify a techrep label if a value other than 1 is desired.

assembly
,,,,,,,,
| Required (no default)
| 
| The name of the genome assembly for the sample, such as ``hg38``.


reference
,,,,,,,,,
| Required (no default, but can be downloaded automatically)
|
| The reference genome file for the sample. Hich will automatically download a genome reference if not provided for the following assemblies:

- ``hg38``, ``homo_sapiens``, ``GRCh38``
- ``mm10``
- ``dm6``
- ``galGal5``
- ``bGalGal5``
- ``danRer11``

chromsizes
,,,,,,,,,,
| Required (no default, but can be built automatically)
|
| The chromsizes file for the `reference`_ genome, a two-column list of contig names and contig sizes in bp. Built automatically from the reference genome if not specified for a given sample.

min_mapq
,,,,,,,,
| Not required
| Default: ``30``
| 
| Minimum alignment threshold to keep an aligned read.

datatype
,,,,,,,,
| Required
| 
| Options:

- ``fastq`` `default`
- ``sam``
- ``bam``
- ``pairs``

| 
| The format for input read data.
| 
| `Note: Hich can read data compressed in gzip format, but gzip compression does not need to be explicitly specified.`

id
,,
| Required (defaults to ``{condition}_{biorep}_{techrep}``)
| 
| A unique id label for the sample.

Alignment
.........

aligner
,,,,,,,
| Required if ``datatype == fastq``
| Available options:

- ``bwa``
- ``bwa-mem2`` `default`

| 
| While `bwa-mem2` is 1.3-3x faster, indexing genomes with ``bwa-mem2`` requires a 60-80 Gb memory footprint, whereas indexing with `bwa` can be done in less than 32 Gb.

index_dir
,,,,,,,,,
| Not required
|
| Directory where the `aligner`_-specific `reference`_ genome index files are stored. Each file should start with the same `index_prefix`_. If not specified, Hich will attempt to index the reference genome and will output the result to `resources/.hich` under a subdirectory for the specific aligner.

index_prefix
,,,,,,,,,,,,
| Not required
|
| Prefix shared by all needed `aligner`_-specifi `reference`_ genome index files in the `index_dir`_ directory. If not specified, Hich will attempt to index the reference genome and will output the result to `resources/.hich` under a subdirectory for the specific aligner.

aligner_threads
,,,,,,,,,,,,,,,
| Default: ``10``
| 
| Max threads to use for alignment. It is highly recommended to set this to the maximum number of available cores. Note that only one alignment process is spawned at a time. This is because every `aligner`_ Hich uses (BWA MEM and BWA MEM2) are internally parallelized, so there is no substantial performance gain to running multiple alignment processes in parallel, while the substantial memory footprint is duplicated for each aligner instance being run.

bwa_flags
,,,,,,,,,
| Default: ``-SP5M``
| 
| Flags to use for the `aligner`_ ``bwa mem`` or ``bwa-mem2 mem``. The default ``-SP5M`` is recommended by `4DN <https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline>`_ for aligning paired-end Hi-C reads with ``bwa mem`` or ``bwa-mem2 mem``. See `bwa manual reference page <https://bio-bwa.sourceforge.net/bwa.shtml>`_ for additional options.

Pairs processing
................

enzymes
,,,,,,,
| Default: none
| 
| If restriction enzymes were used to digest the sample, they can be listed here. Hich allows specifying "Arima" for the Arima Hi-C+ kit enzymes. Any enzymes or combination of enzymes in Biopython's `Bio.Restrict <http://biopython.org/DIST/docs/cookbook/Restriction.html>`_ library can be used. Multiple enzymes should be separated by commas `,`. If specified, a "fragment index" (a digest of the reference genome using the enzymes in .bed format) will be produced automatically, used to tag tne ends of each read with the restriction fragment it maps to, and then filter out any reads where each end maps to the same restriction fragment. If not specified, none of these steps occur. See `fragfile`_ for how to use an already-created fragment index.

.. csv-table:: Example samples.csv specifying different ``enzymes`` options
    :file: tables/reference_enzymes.tsv
    :delim: tab
    :header-rows: 1

fragfile
,,,,,,,,
| Default: none
| 
| An already-created fragment index in .bed format, to be used for tagging contacts with the fragment from which each end originated if the `enzymes`_ parameter is specified for the sample.

deduplicate
,,,,,,,,,,,
| Options:

- true `default`
- false

| 
| Whether to remove technical duplicates (i.e. PCR or optical duplicates). Deduplication is applied to biological replicates after forming them from non-deduplicated technical replicates or after ingesting them directly into Hich. Hich deduplicates technical replicates `after` using them to merge biological replicates.

pairs_format
,,,,,,,,,,,,

chrom1
^^^^^^
| Required
| Default: 2
| 
| The column in the .pairs file where the first chromosome is labeled for each read.

pos1
^^^^
| Required
| Default: 3
| 
| The column in the .pairs file where the first base pair position is labeled for each read.

chrom2
^^^^^^
| Required
| Default: 4
| 
| The column in the .pairs file where the second chromosome is labeled for each read.

pos2
^^^^
| Required
| Default: 5
| 
| The column in the .pairs file where the second base pair position is labeled for each read.

parse_params
,,,,,,,,,,,,
| Default:

- ``--flip``
- ``--drop-readid``
- ``--drop-seq``
- ``--drop-sam``

| 
| Extra parameters to use for parsing .sam/.bam alignments into .pairs format.
| 
| `Note:` The ``drop-*`` parameters are one of the most impactful for making Hich fast and giving it a low disk footprint. It is not recommended to remove these parameters unless you know what you are doing, although additional parameters can be added.

dedup_params
,,,,,,,,,,,,
| Extra parameters to use during the deduplication step.

select_params
,,,,,,,,,,,,,
| Extra parameters to use during the selection step.

select_condition
,,,,,,,,,,,,,,,,
| Read-level filters to use during the selection step.

keep_pair_types
^^^^^^^^^^^^^^^
| Default: ``UU``, ``UR``, ``RU``
| 
| U is for a unique aligned read, whereas an R is "rescued" by detecting pairs where one side maps to locus 1 and the other to a slightly different position on locus 1 and to locus 2, the classic "split ligation junction" pattern that represents an observed, rather than inferred, ligation junction.

keep_trans
^^^^^^^^^^^^^^^
| Options:

- ``true`` `default`
- ``false``

| 
| Whether to keep interchromosomal ("trans") contacts. Note that this should be left as true if forming .mcool files and using the default ``trans-only`` option, which normalizes contact matrices based exclusively on trans contacts, which are in some cases thought to yield more biologically representative results.

keep_cis
^^^^^^^^^^^^^^^
| Options:

- ``true`` `default`
- ``false``

| 
| Whether to keep intrachromosomal ("cis") contacts.

min_dist_fr
^^^^^^^^^^^^^^^
| Default: ``1000``
| Minimum insert size (in bp) to keep FR (+- or inward) strands. In Hi-C, the set of short-range FR strands can be highly enriched in undigested chromatin, which shows up in Hich's MultiQC report as a percentage of FR orientations substantially higher than the expected 25%. These can be filtered out using this option.

min_dist_rf
^^^^^^^^^^^^^^^
| Default: ``1000``
| Minimum insert size (in bp) to keep RF (-+ or outward) strands. In Hi-C, the set of short-range FR strands can be highly enriched in self-circles (digested fragments that self-ligated end to end), which shows up in Hich's MultiQC report as a percentage of RF orientations substantially higher than the expected 25%. These can be filtered out using this option.

min_dist_ff
^^^^^^^^^^^^^^^
| Default: ``0``
| Minimum insert size (in bp) to keep FF (++) strands.

min_dist_ff
^^^^^^^^^^^^^^^
| Default: ``0``
| Minimum insert size (in bp) to keep RR (--) strands.

chroms
^^^^^^^^^^^^^^^
| If specified, each read alignment must be to a chromosome in this set.

discard_same_frag
^^^^^^^^^^^^^^^^^^
| Options:

- ``true`` `default`
- ``false``

| 
| If true, fragments whose alignments are mapped to restriction fragments will be discarded if both ends mapped to the same restriction fragment.

Matrix processing
.................

make_hic
,,,,,,,,,,,
| Arguments supplied to juicer tools' ``pre`` command when forming a Hi-C contact matrix.

make_cool
,,,,,,,,,,
| Arguments supplied to the ``cooler cload`` command for forming .cool format precursors to the .mcool contact matrix.

make_mcool
,,,,,,,,,,
| Default:

- ``--balance``
- ``--balance-args 'max-iters 2000 --trans-only'``

| 
| Arguments supplied to the ``cooler zoomify`` command for coarsening high-res .cool matrices into multi-resolution .mcool contact matrices. The chosen defaults will generate multi-res contact matrices containing both the raw contacts and balancing weights produced using the trans contacts only.

matrix
,,,,,,

make_mcool_file_format
^^^^^^^^^^^^^^^^^^^^^^
| Options:

- true `default`
- false

| 
| Whether to produce .mcool-format contact matrices (the Open2C multi-resolution format). Currently required for feature calling and QC.

make_hic_file_format
^^^^^^^^^^^^^^^^^^^^^^
| Options:

- true `default`
- false `default`

| 
| Whether to produce .hic-format contact matrices (compatible with the Juicer tool ecosystem including the Juicebox browser).

resolutions
^^^^^^^^^^^^^^^^^^^^^^
| Default:

- 1000
- 2000
- 5000
- 10000
- 20000
- 50000
- 100000
- 200000
- 500000
- 1000000

| 
| Reference chromosome coordinates will be partitioned into these uniform block sizes (in bp) and contact ends mapped to those blocks to generate contact matrices. Lower numbers represent higher-resolution matrices. 

Quality control
...............

hicrep
,,,,,,

call_on
^^^^^^^
| Options:

- is_techrep `default`
- is_biorep `default`
- is_condition `default`

| 
| Whether to compute Hicrep SCC scores on technical replicates, biological replicates, and conditions. Results for all comparisons are output to a single .tsv file with a per-column header giving the pair of samples, chromosome, resolution, and Hicrep parameters that were used, along with the SCC score.

resolutions
^^^^^^^^^^^
| Default:

- 10000
- 100000
- 1000000

| 
| Which resolutions to use for calling Hicrep SCC scores.

chroms
^^^^^^
| Which chromosomes to use for calling Hicrep SCC scores. If not specified, all chromosomes shared by both matrices at the given resolution will be used.

exclude
^^^^^^^
| Which chromosomes to exclude when calling Hicrep SCC scores.

chrom_filter
^^^^^^^^^^^^
| A conditional statement in Python to determine whether to use a chromosome for Hicrep as a function of its name (referenced via the ``chrom`` variable) and size (the ``size`` variable). It will be evaluated using Python's ``eval`` statement.

h
^^^^^
| Values of Hicrep's ``h`` parameter to use.

dBPMax
^^^^^^
| Values of Hicrep's ``dBPMax`` parameter to use.

bDownSample
^^^^^^^^^^^^
| Values of Hicrep's ``bDownSample`` parameter to use.

Feature calling
...............

compartments
,,,,,,,,,,,,

resolution
^^^^^^^^^^^
| Default: 5000
| 
| The resolution at which compartments should be called.

cooltools_eigs_cis_params
^^^^^^^^^^^^^^^^^^^^^^^^^^
| Defaults:

- --bigwig

| Additional parameters that should be passed to cooltools_eigs_cis. The default specifies that a .bigwig-format file should be generated as well as the .bedgraph format.

insulation
,,,,,,,,,,

resolution
^^^^^^^^^^^
| Default: 5000
| 
| The resolution at which insulation should be called.

cooltools_insulation_params
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| Defaults:

- --bigwig

| Additional parameters that should be passed to cooltools_eigs_cis. The default specifies that a .bigwig-format file should be generated as well as the .bedgraph format.

loops
,,,,,,,
| Hich uses Mustache for loop and differential loop calling. This software was chosen mainly for its theoretical advantages. Based on scale space theory, it applies an artifact-free filter to a matrix to remove fine details, then detects blobs which are called as loops. It thereby takes local information into account in loop calling. Differential loops for a pair of matrices are loops that are present or enriched in one matrix and not present or depleted in the other. An added practical benefit is that Mustache is fast enough that it can run on a CPU, whereas many other loop callers require a GPU.

call_on
^^^^^^^^^^^
| Options:

- is_techrep `default`
- is_biorep `default`
- is_condition `default`

| 
| Whether to call loops on technical replicates, biological replicates, and conditions.

use_format
^^^^^^^^^^^^^^
| Options:

- mcool `default`
- hic

| Mustache can use both .mcool and .hic matrix formats as input. Loops will only be called on samples where the appropriate matrix type is output. If both are generated, which is used should not affect the outcome.

mustache_params
^^^^^^^^^^^^^^^^^
| Default:

- ``--resolution 5000``
- ``--pThreshold .1``
- ``--sparsityThreshold .88``

| 
| Parameters passed to mustache_diffloops, which will output both individual matrix loop calls and a pair of diffloops calls for each matrix.