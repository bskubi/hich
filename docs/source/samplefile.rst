Sample File
===========

The Hich sample file is a tab- or comma-delimited (.tsv or .csv) file that describes your experimental setup. It can often look as simple as this:

Example 1. Multiple conditions, bioreps, techreps, and reference genomes
----------------------------------------------------------------------------

.. csv-table::
    :file: tables/reference_samplefile1.tsv
    :delim: tab
    :header-rows: 1

The condition, biorep, and techrep columns label the condition, biological replicate, and technical replicate for each sample (row). The fastq1 and fastq2 columns give the path to the paired-end fastq files containing sequencing data for each sample. The assembly column contains the name of the reference genome (hg38 = human, mm10 = murine). The enzymes column can be a comma-separated list of restriction enzymes or the name of the kit used in the experiment.

With this setup, Hich will automatically download the reference genome (provided it recognizes the assembly name), indexing it for alignment with bwa-mem2, and producing other needed resource files such as a partition of the genome into restriction fragments. Note that resource files will be produced only once (so that here, one copy of the hg38 reference will be downloaded and shared by all samples).

Example 2. A different aligner
---------------------------------

.. csv-table::
    :file: tables/reference_samplefile4.tsv
    :delim: tab
    :header-rows: 1

Here, only one sample is being used. The biorep and techrep will default to "1". There is also a new column, specifying that the "bwa" aligner should be used rather than the default "bwa-mem2". While the default bwa-mem2 aligner is faster, the alternative bwa aligner requires less memory for indexing, which is the main motivation for offering it as an option. Both aligners should produce identical outputs.

Example 3. Use local reference, index, chromsizes, and fragfile
---------------------------------------------------------------------------------------------------
.. csv-table::
    :file: tables/reference_samplefile5.tsv
    :delim: tab
    :header-rows: 1

Here, all resource files are already created and their paths are given. The new columns are:

- reference: The reference genome in FASTA format
- chromsizes: A two-column tab-delimited file with the contigs in column 1 and the size (bp) of the contig in column 2.
- index_dir: The directory where the desired index files are located
- index_prefix: The common prefix string for each index file. BWA, for example, uses the files [index_prefix].amb, [index_prefix].ann, [index_prefix].bwt, [index_prefix].pac, and [index_prefix].sa.
- fragfile: A partition of the genome into restriction fragments in BED format.

The ``index_dir`` and ``index_prefix`` columns must either both have a value or neither have a value, but other than that, it is fine to specify some resource file columns and leave others blank for any given sample.

4D Nucleome makes several commonly-used BWA index files available for manual download. If using them, you will need to download and unzip them before putting their paths into the sample file. Many labs produce an index for the genomes they work with and reuse them across many projects. As indexing the genome is the most RAM-intensive process Hich executes, it can be helpful to supply one if available.

Links to 4DN BWA Index files:
- `danRer11 (compressed) <https://data.4dnucleome.org/files-reference/4DNFIUH46PG1/#details>`_
- `mm10 (uncompressed) <https://data.4dnucleome.org/files-reference/4DNFIZ2PWCC2/#details>`_
- `hg38 (uncompressed) <https://data.4dnucleome.org/files-reference/4DNFIZQB369V/#details>`_
- `galGal5 (compressed) <https://data.4dnucleome.org/files-reference/4DNFIVGRYVQF/#details>`_
- `dm6 (compressed) <https://data.4dnucleome.org/files-reference/4DNFIO5MGY32/#details>`_
- `mm10 (compressed) <https://data.4dnucleome.org/files-reference/4DNFI823LSI8/#details>`_

If using one of these, it is essential to set the ``aligner`` column to ``bwa`` in the sample file for each sample using a BWA index, since Hich defaults to using the faster (but more RAM-intensive) bwa-mem2, which uses a different index format.

Merging and downsampling
------------------------

In this documentation, some terms have specific definitions.

- **Sample** A sample is the complete collection of data and parameters associated with a single technical replicate (techrep), biological replicate (biorep), or condition.
- **Biological sample (biosample)** A biosample refers to the physical material (cells, a sequencing library prep) that is associated with a particular Hich sample.
- **Condition** An experimental condition, such as "No treatment (NT)" or "Knockout (KO)" or data obtained from different cell lines or drug exposures. Conditions are formed by merging their respective bioreps, but are not deduplicated as their input biorep samples will have already been deduplicated.
- **Biological replicate** Each condition can have one or more biological replicates, representing data from the same experimental condition applied to one or more biosamples. Bioreps are formed by merging their respective techreps, and are then deduplicated before contact matrices are generated. 
- **Technical replicate** Each biorep can have one or more technical replicates, representing data from the same biorep sequenced multiple times. After merging techreps to form bioreps (which are independently deduplicated), the individual techreps are deduplicated as well. 
- **Input sample** A sample that is input to Hich as a row in the sample file
- **Synthetic sample** A sample that is created by Hich based on input samples, typically via merging and/or downsampling.

When Hich merges techrep samples to form biorep samples, it keeps the techrep samples as well. It also keeps the biorep samples after merging them to form condition samples. That means that the output contains techrep, biorep, and condition samples for the experiment, facilitating downstream QC and analysis.

This is a test.