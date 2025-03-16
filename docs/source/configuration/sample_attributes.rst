Sample attributes reference
...........................


Relationships between samples
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

``id``
    String, required, must be different for each sample. Used as filename prefix for output files.

.. note::
    ``id`` is typically built algorithmically by concatenating ``techrep``, ``biorep``, and ``condition`, as well as the aggregate profile name, but can be manually specified.

``condition``
    Optional. Labels basal tier in the experimental design hierarchy.

``biorep``
    Optional. Labels secondary tier in experimental design hierarchy. Samples with distinct ``techrep`` and ``biorep`` but common ``condition`` can be aggregated.

``techrep``
    Optional. Labels tertiary tier in experimental design hierarchy. Samples with distinct ``techrep`` but common ``condition`` and ``biorep`` can be aggregated.

Resource files
,,,,,,,,,,,,,,

``assembly``
    Genome assembly label.

``genomeReference``
    Path or URL to genome reference fasta file. If ``genomeReference`` unspecified but ``assembly`` is one of the supported options, Hich downloads the genome reference from the ENCODE project or NCBI. If multiple samples will use the downloaded reference, it will only be downloaded once and shared by all the samples that need it. Supported options for automatically downloading genome reference: ``hg38``, ``homo_sapiens``, or ``GRCh38`` (GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz), ``mm10`` (mm10_no_alt_analysis_set_ENCODE.fasta.gz), ``dm6`` (GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz), ``galGal5`` or ``bgalGal5`` (GCA_027408225.1_bGalGal5.pri_genomic.fna.gz), ``danRer11`` (GCF_000002035.6_GRCz11_genomic.fna.gz), or ``ce10`` (GCF_000002985.6_WBcel235_genomic.fna.gz).

``indexDir``
    Path to directory containing aligner index.

``indexPrefix``
    Common prefix of all aligner index files (it is required that all aligner index files share a common prefix). Example: if the aligner index files are ``hg38.*``, then ``indexPrefix`` should be ``hg38``.

``chromsizes``
    Path to tab-delimited headerless file with contig names in first column, length of contig in base pairs as second column. Automatically created based on ``genomeReference`` if unspecified and shared among samples with a common reference that all left ``genomeReference`` unspecified.

``restrictionEnzymes``
    Optional. Space-delimited list of restriction enzyme names used in restriction digest for the sample.

``fragmentIndex``
    Optional. Path to BED file containing start and end positions of restriction fragments for the digest used for the sample. If the ``restrictionEnzymes`` option is specified but ``fragmentIndex`` is not, then Hich will create a ``fragmentIndex`` file based on the ``restrictionEnzymes`` and ``genomeReference`` and share it among samples with the same reference and enzymes. Any combination of enzymes in the REBASE database as accessed via `biopython's restriction enzymes module <https://biopython.org/DIST/docs/cookbook/Restriction.html>`_ can be used, as well as ``Arima``, ``Phase Proximo 2021+ Plant`` (or ``Phase Plant``), ``Phase Proximo 2021+ Animal`` (or ``Phase Animal``), ``Phase Proximo 2021+ Microbiome`` (or ``Phase Microbiome``), ``Phase Proximo 2021+ Human`` (or ``Phase Human``) or ``Phase Proximo 2021+ Fungal`` (or ``Phase Fungal``). 

Aligning reads
,,,,,,,,,,,,,,,,,,,,

Hich toolkit: `bwa mem <https://github.com/lh3/bwa>`_, `bwa-mem2 <https://github.com/bwa-mem2/bwa-mem2>`_, `bsbolt <https://bsbolt.readthedocs.io/en/latest/>`_

See also: :ref:`Resource files` under ``assembly``, ``genomeReference``, ``indexDir``, ``indexPrefix``.

``fastq1``
    Path to R1 or single-end read fastq file.

``fastq2``
    Path to R2 fastq file. Leave blank or unspecified if using single-end reads.

``aligner``
    Aligner to use for aligning the sample. Options: ``bwa`` (slower, lower memory footprint), ``bwa-mem2`` (fast, higher memory footprint), ``bsbolt`` (methyl Hi-C)

``bwaFlags``
    CLI options passed to aligner (note that all aligners including BSBolt are based on ``bwa mem``). Typically, use ``-SP5M``. Do not use the ``bwa`` option ``-t`` or the BSBolt options ``-OT``, ``-O``, ``-DB``, ``-F1``, or ``-F2`` as these are hardcoded by Hich based on other sample attributes.

``minMapq``
    Reads below this MAPQ cutoff will be discarded. Note that different aligners approximate MAPQ differently. The approach used by ``bwa`` is what's relevant for Hich.

Hi-C contacts ingested from .pairs or parsed from .sam/.bam
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

Hich toolkit: `samtools <https://github.com/samtools/samtools>`_, `pairtools <https://pairtools.readthedocs.io/en/latest/cli_tools.html>`_

``sambam``
    Path to .sam/.bam file containing aligned reads which will be parsed using ``pairtools parse2`` to obtain a `4DN .pairs <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_ file. 

``pairs``
    Path to `4DN .pairs <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_ file containing Hi-C contacts to ingest.

``latestPairs``
    Path to `4DN .pairs <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_ file containing Hi-C contacts to ingest or the result of most recent data processing step on the sample contacts.

``pairtoolsParse2Params``
    List of flags passed to `pairtools parse2 <https://pairtools.readthedocs.io/en/latest/cli_tools.html#pairtools-parse2>`_. Uses ``minMapq`` if specified for the sample as a default value for the ``pairtools parse2 --min-mapq`` option, but this can be overridden by passing ``--min-mapq`` n ``pairtoolsParse2Params``. Hardcoded options that should not be provided here: ``--flip``, ``--assembly``, ``--chroms-path``.

.. note::
    Hich will inspect .sam/.bam files to determine if they are sorted, and sort them automatically by name (required for inputs to ``pairtools parse``) only if necessary. It will then sort the output by position.

.. note::
    ``pairtools parse2`` has a ``--drop-readid`` parameter, which can drastically shrink the disk space required for the .pairs file. This is useful, but for single cell data (see below), it was challenging to engineer a way to drop this column when it's necessary to extract the ``cellID`` column value from the ``readID`` column of the .sam/.bam file used as input to parsing the .pairs file. For this reason, the ``--drop-readid`` parameter is not actually passed to ``pairtools parse2``. Instead, ``--placeholder readID .`` is passed to ``hich reshape``, which accomplishes the same result while permitting ``cellID`` to be extracted from the ``readID`` column if necessary.

Optional single-cell attributes
_______________________________

.. note::
    These attributes can be ignored for bulk data. For single cell-aware fragment filtering, deduplication and to maintain cell ID for future analysis, Hich must put a unique identifier for the cell attributed to each contact in the .pairs file into a column labeled ``cellID`` in the .pairs file. This identifier can be extracted by Hich automatically from the read ID or from a .sam/.bam tag using the sample attributes in this section using the Hich CLI command ``hich reshape``.

``cellBarcodeField``
    Required if parsing cell ID from a .sam/.bam file. Should be either ``readID`` the name of a .sam/.bam tag. This field will be parsed for each read in the .sam/.bam file in order to extract the value of the ``cellID``. The patterns used to accomplish this extraction are specified below.

``cellBarcodeRegexPattern``
    Optional. Should be a Python regex compatible with re (regexes can be tested at `regex101.com <https://regex101.com/>`_). Along with ``cellBarcodeGroup``, the regex will be applied to parse the field specified in ``cellBarcodeField`` and the match will be put into the ``cellID`` field of the .pairs file. Overrides ``cellBarcodeParsePattern`` if both are specified.

``cellBarcodeGroup``
    Optional. An integer specifying which match group from the regex specified by ``cellBarcodeRegexPattern`` should be used as the value of ``cellID``. 0 uses all match groups. Defaults to 0 if ``cellBarcodeField`` and ``cellBarcodeRegexPattern`` are specified and ``cellBarcodeGroup`` is not.

``cellBarcodeParsePattern``
    Optional. An alternative and potentially simpler way to parse ``cellBarcodeField`` by using Python's `parse <https://pypi.org/project/parse/>`_ library syntax. From the pattern specified the ``{cellID}`` named part will be extracted and put into the ``cellID`` column in the .pairs file. Example: ``{}:{cellID}`` will extract the part after a colon (:) and put it into the ``cellID`` column.

``globalDefaultReshapeToCellID``
    Optional. Must be specified in the params file or nextflow.config. If ``cellBarcodeField`` is specified for a sample but either ``cellBarcodeRegexPattern`` nor ``cellBarcodeParsePattern`` is specified, then ``globalDefaultReshapeToCellID`` is used to determine how the ``cellID`` column will be parsed. Ignored if ``cellBarcodeRegexPattern`` or ``cellBarcodeParsePattern`` is given for the sample. 

``globalDefaultReshapeToCellID.option``
    Optional. Either ``--regex`` or ``--parse``, which determines whether ``globalDefaultReshapeToCellID.pattern`` (below) will be parsed using Python's ``re`` library or its ``parse`` library (see above options for details).

``globalDefaultReshapeToCellID.pattern``
    Optional. Interpreted either a Python ``re`` regex or Python ``parse`` pattern depending on the value of ``globalDefaultReshapeToCellID.option``.

``globalDefaultReshapeToCellID.group``
    Optional. The match group to use for the regex. Ignored if unspecified, and should be left unspecified if using ``parse``.

``reshapeParams``
    Optional additional params passed to ``hich reshape``.

Filtering Hi-C contacts
,,,,,,,,,,,,,,,,,,,,,,,

Hich toolkit: `pairtools <https://pairtools.readthedocs.io/en/latest/cli_tools.html>`_

See also: :ref:`Resource files` under ``restrictionEnzymes``, ``fragmentIndex``

``selectFilters``
    A multi-attribute of filters to apply to Hi-C contacts in .pairs files.

``selectFilters.keepPairTypes``
    Pairtools `pair types <https://pairtools.readthedocs.io/en/latest/parsing.html>`_ to keep. Keeping ``UU``, ``UR``, and ``RU`` is recommended.

``selectFilters.keepTrans``
    If false, discards reads mapping to different chromosomes/contigs. If unspecified, these contacts will be kept.

``selectFilters.keepCis``
    If false, discards reads mapping to the same chromosome/contig. If unspecified, these contacts will be kept.

``selectFilters.minDistFR``
    If specified, then for reads with the orientation FR, discards if they are below this distance between ``pos1`` and ``pos2``.

``selectFilters.minDistRF``
    If specified, then for reads with the orientation RF, discards if they are below this distance between ``pos1`` and ``pos2``.

``selectFilters.minDistFF``
    If specified, then for reads with the orientation FF, discards if they are below this distance between ``pos1`` and ``pos2``.

``selectFilters.minDistRR``
    If specified, then for reads with the orientation RR, discards if they are below this distance between ``pos1`` and ``pos2``.

.. note::
    Two technical artifacts that routinely appear in Hi-C experiments enriched in short-range contacts are undigested chromatin and self-ligated strands. These will appear in the multiQC reports generated by Hich as a strong enrichment in the FR and RF orientations below a certain distance threshold. By pausing the Hich run after parsing to pairs and inspecting this report, the ``minDist`` values can be chosen appropriately according to the QC data. Data with no strand bias should have very close to 25% of each orientation.

``selectFilters.discardSingleFrag``
    Discard contacts where both ends map to the same restriction fragment as these likely originate from undigested chromatin. Requires that samples have been tagged with this information, which Hich will do automatically if ``fragmentIndex`` is specified.

``pairtoolsSelectParams``
    Additional parameters to pass to ``pairtools select``. The following options are hardcoded in Hich and should not be specified here: ``--output-rest``, ``--output``, ``--nproc-in``, ``--nproc-out``.

Downsampling, merging, and deduplicating samples
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

Creating contact matrices
,,,,,,,,,,,,,,,,,,,,,,,,,

Generating multiQC reports on Hi-C contacts
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

Calling HiCRep SCC scores
,,,,,,,,,,,,,,,,,,,,,,,,,

Calling compartment scores
,,,,,,,,,,,,,,,,,,,,,,,,,,

Calling insulation scores
,,,,,,,,,,,,,,,,,,,,,,,,,

Calling loops
,,,,,,,,,,,,,,,,

Calling differential loop enrichments (diffloops)
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

Recent outputs
,,,,,,,,,,,,,,

``latest``

``latestSambam``

``latestPairs``

``latestMatrix``

Hich sample attributes used internally (not typically specified by user)
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

``isSingleCell``