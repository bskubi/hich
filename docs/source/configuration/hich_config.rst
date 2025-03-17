Hich pipeline options reference
...................................

Typically specified in :ref:`params file <Params file>`
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

``defaults``
    Hashmap. Required. :ref:`Sample attributes <Sample attributes reference>` specified in the ``defaults`` block are applied to all samples by default, unless overridden on for individual samples, (i.e. in the sample file). YAML example:

.. code:: c

    defaults:
        minMapq: 30
        aligner: "bwa-mem2"
        bwaFlags: ["-SP5M"]
        pairtoolsSelectFilters:
            keepPairTypes: ["UU", "UR", "RU"]
            keepTrans: true
            keepCis: true

``aggregate``
    Hashmap. Optional. :ref:`Aggregate profiles <Downsampling, merging, and deduplicating samples>` are declared here and used control whether and how samples will be downsampled, merged and deduplicated. The keys of the aggregate hashmap are profile names. Values are hashmaps defining how that profile behaves. YAML example:

.. code:: c

    aggregate:
        profile1:
            dedupMaxMismatch: 3
            dedupMethod: "max"
            techrepDedup: true
        profile2:
            mergeTechrepToBiorep: true

``sampleSelectionStrategies``
    Hashmap. Optional. :ref:`Sample selection strategies` are declared here and used to control the samples that are used as inputs into various feature calling methods. The keys of the ``sampleSelectionStrategies`` hashmap are strategy names. Values are a hashmap of sample attributes and lists of acceptable values of those attributes. YAML example:

.. code:: c

    sampleSelectionStrategies:
        strategy1:
            condition: ["c1", "c2"]
            aggregateProfile
        strategy2:

``hicrep``
    Hashmap. Optional. :ref:`Calling HiCRep SCC scores` is controlled here. Keys are names of parameterization profiles for hicrep. Values are analysis plans, which are the parameterizations themselves.

``loops``
    Hashmap. Optional. :ref:`Calling loops` is controlled here. Keys are names of parameterization profiles for Mustache. Values are analysis plans, which are the parameterizations themselves.

``differentialLoops``
    Hashmap. Optional. :ref:`Calling differential loop enrichments (diffloops)` is controlled here. Keys are names of parameterization profiles for Mustache. Values are analysis plans, which are the parameterizations themselves.

``compartments``
    Hashmap. Optional. :ref:`Calling compartment scores` is controlled here. Keys are names of parameterization profiles for Mustache. Values are analysis plans, which are the parameterizations themselves.

``insulation``
    Hashmap. Optional. :ref:`Calling insulation scores` is controlled here. Keys are names of parameterization profiles for Mustache. Values are analysis plans, which are the parameterizations themselves.


Typically specified in :ref:`nextflow.config <Nextflow config file>`
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

``sampleFileSep``
    Single-character string. Required to parse sample file. Column separator. Use "\t" for tab (TSV files) or "," for comma (CSV files). Other settings can be used as well.

``humid``
    Boolean. Required. If true, then ingested, gzipped fastq files will be downsampled to the number of reads specified in ``general.humidDefault``.

``general``
    Hashmap. Required. Contains additional parameters.

``general.humidDefault``
    List of strings. Required. Specified in provided nextflow.config. For gzipped fastq files, the number of reads to downsample to for a "humid" run.

``general.hichContainer``
    String. Required. Specified in provided nextflow.config. Location of Hich CLI utilities Docker container.

``general.chromsizesContainer``
    String. Required. Specified in provided nextflow.config. Location of ucsc-fasize Docker container used to produce chromsizes file from genome reference.

``general.mustacheContainer``
    String. Required. Specified in provided nextflow.config. Location of Docker container used to call Mustache loops and differential loops.

``general.juicerContainer``
    String. Required. Specified in provided nextflow.config. Location of Docker container used to call Juicer tools to produce .hic contact matrix.

``general.hictkContainer``
    String. Required. Specified in provided nextflow.config. Location of Docker container used to convert between .mcool and .hic formats.

``general.qcAfter``
    List of strings. Required. Specified in provided nextflow.config. Processes after which multiQC reports should be generated.

``general.publish``
    Hashmap. Required. Contains additional parameters used to control where outputs are published.

``general.publish.mode``
    String. Required. Mode used to publish outputs. See `Nextflow publishDir documentation for options <https://www.nextflow.io/docs/latest/reference/process.html#publishdir>`_

``general.publish.genomeReference``
    String. Required if downloading genome reference. Target directory where downloaded genome references will be published.

``general.publish.chromsizes``
    String. Required if auto-producing chromsizes. Target directory where chromsizes file will be published.

``general.publish.bwaMem2Index``
    String. Required if auto-producing bwa-mem2 aligner index. Target directory where bwa-mem2 aligner index will be published.

``general.publish.bwaIndex``
    String. Required if auto-producing bwa aligner index. Target directory where bwa aligner index will be published.

``general.publish.bsboltIndex``
    String. Required if auto-producing bsbolt aligner index. Target directory where bsbolt aligner index will be published.

``general.publish.fragmentIndex``
    String. Required if auto-producing restriction digest fragment index. Target directory where fragment index will be published.

``general.publish.align``
    String. Required if aligning .fastq files. Target directory where sam/bam files will be published.

``general.publish.parse``
    String. Required if parsing sam/bam files to .pairs format. Target directory where resulting .pairs files will be published.

``general.publish.dedup``
    String. Required if deduplicating .pairs format files. Target directory where resulting .pairs files will be published.

``general.publish.mcool``
    String. Required if generating .mcool files either from .pairs files or by converting from .hic format. Target directory where resulting .mcool files will be published.

``general.publish.hic``
    String. Required if generating .hic files either from .pairs files or by converting from .mcool format. Target directory where resulting .hic files will be published.

``general.publish.pairStats``
    String. Required if generating pairtools stats files. Target directory where resulting stats files will be published.

``general.publish.qc``
    String. Required if generating multiQC reports. Target directory where resulting reports will be published.