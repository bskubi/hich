Hich Reference
==============

In Hich, a **sample** is a single unit of data, such as a **technical replicate**, **biological replicate**, or experimental **condition**.

Sample attributes
-----------------

General
.......

sample_id
,,,,,,,,,
Required (no default)

A unique name for the sample.

id
,,
| Required (defaults to sample_id)
| 
| A unique name for the sample. Automatically chosen by Hich, not set manually.

assembly
,,,,,,,,
| Required (no default)
| 
| The genome assembly for the sample, such as ``hg38``. Hich will automatically download a genome reference if not provided for the following assemblies:

- ``hg38``, ``homo_sapiens``, ``GRCh38``
- ``mm10``
- ``dm6``
- ``galGal5``
- ``bGalGal5``
- ``danRer11``

datatype
,,,,,,,,
| Required (no default)
| 
| Options:

- ``fastq``
- ``sam``
- ``bam``
- ``pairs``

| 
| The format for input read data.
| 
| `Note: Hich can read data compressed in gzip format, but gzip compression does not need to be explicitly specified.`


is_techrep
,,,,,,,,,,
| Required
| Available options:


- ``true`` `default`
- ``false``

| 
| Whether or not the sample is a technical replicate.

is_biorep
,,,,,,,,,
| Required
| Available options:

- ``true``
- ``false`` `default`

| 
| Whether or not the sample is a biological replicate.

is_condition
,,,,,,,,,,,,
| Required
| Available options:

- ``true``
- ``false`` `default`

| 
| Whether or not the sample is a condition.

min_mapq
,,,,,,,,
| Not required
| Default: ``30``
| 
| Minimum alignment threshold to keep an aligned read.

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

aligner_threads
,,,,,,,,,,,,,,,
| Default: ``10``
| 
| Max threads to use for alignment. It is highly recommended to set this to the maximum number of available cores. Note that only one alignment process is spawned at a time. This is because the aligners Hich uses (BWA MEM and BWA MEM2) are internally parallelized, so there is no substantial performance gain to running multiple alignment processes in parallel, while the substantial memory footprint is duplicated for each aligner instance being run.

bwa_flags
,,,,,,,,,
| Default: ``-SP5M``
| 
| Flags to use for ``bwa mem`` or ``bwa-mem2 mem``. The default ``-SP5M`` is recommended by `4DN <https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline>`_ for aligning paired-end Hi-C reads with ``bwa mem`` or ``bwa-mem2 mem``. See `bwa manual reference page <https://bio-bwa.sourceforge.net/bwa.shtml>`_ for additional options.

Pairs processing
................

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