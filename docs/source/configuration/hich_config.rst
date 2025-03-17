Hich-config reference
...................................

.. _params-file:

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
    Hashmap. Optional. :ref:`Sample selection strategies` are declared here and used to control the samples that are used as inputs into various feature calling methods. The keys of the ``sampleSelectionStrategies`` hashmap are strategy names. Values are a hashmap of sample attributes and lists of acceptable values of those attributes, which can be given as a single acceptable value or as a list of acceptable values. YAML example:

.. code:: c

    sampleSelectionStrategies:
        strategy1:
            condition: ["c1", "c2"]
            biorep: "b1"
            aggregateProfile: "profile1"
        strategy2:
            condition: "c1"
            aggregateProfile: "profile2"

Analysis methods and analysis plans
-----------------------------------

An analysis plan is a way of parameterizing an analysis method, such as calling HiCRep SCC scores, loops, or compartment scores. They consist of parameter values as well as a sample selection strategy that determines the samples on which the analysis plan will be run. For an analysis plan, its sample selection strategy may be left out, or specified as a single strategy or list of them. If it is a list, the sample must satisfy all the strategies to be included in the analysis plan. If no sample selection strategy is given, all samples will be used as inputs to that analysis plan. Hich manages conflicts by replacing earlier settings with those from more recent strategies, so ``sampleSelectionStrategy: ["strategy1", "strategy2"]`` replaces any attribute requirements in strategy1 that also appear in strategy2 with the required values in strategy2. 

``hicrep``
    Hashmap. Optional. Analysis plans for :ref:`calling HiCRep SCC scores <Calling HiCRep SCC scores>` are specified in this block. Keys are names of analysis plans for hicrep. Values are analysis plans themselves. YAML example:

.. code:: c

    hicrep:
        sccPlan1:
            sampleSelectionStrategy: "strategy1"
            h: 1
            resolutions: [2000, 10000]
            bDownSample: true
        sccPlan2:
            sampleSelectionStrategy: ["strategy1", "strategy2"]
            h: 2
            resolutions: 20000

This defines two analysis plans. The first uses samples conforming to a sample selection strategy named "strategy1". Using the example provided in the reference for sample selection strategies above, samples are required to have a condition label that is either "c1" or "c2", a biological replicate label "b1", the aggregateProfile label "profile1." The second requires that both "strategy1" and "strategy2" match the samples used as input. Using the example, only samples with the condition "c1", biological replicate label "b1", and aggregateProfile label "profile1" will be chosen. All combinations of the other parameters will be used, so ``sccPlan1`` will run on ``(h: 1, resolution: 2000, bDownSample: true)`` and ``(h: 1, resolution: 10000, bDownSample: true)``.

``loops``
    Hashmap. Optional. Analysis plans for :ref:`calling loops <Calling loops>` are specified in this block. Keys are names of analysis plans for hicrep. Values are analysis plans themselves. YAML example:

.. code:: c

    loops:
        loopPlan1:
            sampleSelectionStrategy: "strategy1"
            mustacheParams: ["-r 2000", "-ch chr1"]

Analysis plans simply pass arguments in mustacheParams directly to `Mustache <https://github.com/ay-lab/mustache/blob/master/README.md>`_, so the parameters specified there should be used for these analysis plans. Do not use the ``-f`` or ``-o`` arguments as these are hardcoded into Hich.

``differentialLoops``
    Hashmap. Optional. Analysis plans for :ref:`calling differential loop enrichments (diffloops) <Calling differential loop enrichments (diffloops)>` are specified in this block. Keys are names of analysis plans for hicrep. Values are analysis plans themselves. YAML example:

.. code:: c

    differentialLoops:
        diffloopPlan1:
            mustacheParams: ["-r 2000"]

Analysis plans simply pass arguments in mustacheParams directly to `Mustache diffloops <https://github.com/ay-lab/mustache/blob/master/README.md>`_, so the parameters specified there should be used for these analysis plans. Do not use the ``-f1``, ``-f2``, or ``-o`` arguments as these are hardcoded into Hich.

``compartments``
    Hashmap. Optional. Analysis plans for :ref:`calling compartment scores <Calling compartment scores>` are specified in this block. Keys are names of analysis plans for hicrep. Values are analysis plans themselves. The ``resolution`` value must be specified. YAML example:

.. code:: c

    compartments:
        compartmentPlan1:
            resolution: 2000
            hichCompartmentsParams: ["--chroms chr1,chr2,chr3"]

Analysis plans simply pass arguments in hichCompartmentsParams directly to the Hich CLI utilities ``hich compartments`` method, so the parameters specified there should be used for these analysis plans. Do not pass ``--n_eigs`` as this is hardcoded into Hich (compartment scores based on the first 3 eigenvalues will be generated).

``insulation``
    Hashmap. Optional. Analysis plans for :ref:`calling insulation scores <Calling insulation scores>` are specified in this block. Keys are names of analysis plans for hicrep. Values are analysis plans themselves. YAML example:

.. code:: c

    insulation:
        insulationPlan1:
            resolution: 2000
            cooltoolsInsulationParams: ["--threshold 1"]
            window: 1000

Analysis plans simply pass arguments in cooltoolsInsulationParams will be passed directly to `cooltools insulation <https://cooltools.readthedocs.io/en/latest/cli.html#cooltools-insulation>`_, so the parameters specified there should be used for these analysis plans. Do not pass ``--output`` as this is hardcoded into Hich.

.. _nextflow-config-file:

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