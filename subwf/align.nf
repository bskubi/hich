include {transpack} from './extraops.nf'

process BwaMem2Align {
    publishDir params.general.publish.bam ? params.general.publish.bam : "results",
               saveAs: {params.general.publish.bam ? it : null},
               mode: params.general.publish.mode

    container "bskubi/hich:latest"
    maxRetries 6
    memory {20.GB + 10.GB * (task.attempt-1)}

    // NOTE: Alignment speed is trivially parallelizeable and does not benefit
    // from running alignment in parallel multiple files at once. Each instance
    // of bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
    // tell nextflow to run one alignment process at a time with maxForks 1.

    // NOTE 2: I need to decide on how to work with the possibility that users
    // will run multiple aligners (bwa-mem2, bwa, bsbolt, etc). We should still
    // only run one at a time.
    maxForks 1

    input:
    tuple val(sample_id), path(index_dir), val(index_prefix), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    shell:
    align = "bwa-mem2 mem -t 10 -SP5M ${index_dir}/${index_prefix} ${fastq1} ${fastq2}"
    tobam = "samtools view -b -o ${sample_id}.bam"
    "${align} | ${tobam}"

    stub:
    "touch ${sample_id}.bam"
}

process BwaMemAlign {
    publishDir params.general.publish.bam ? params.general.publish.bam : "results",
               saveAs: {params.general.publish.bam ? it : null},
               mode: params.general.publish.mode

    container "bskubi/hich:latest"
    maxRetries 6
    memory {20.GB + 10.GB * (task.attempt-1)}

    // NOTE: Alignment speed is trivially parallelizeable and does not benefit
    // from running alignment in parallel multiple files at once. Each instance
    // of bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
    // tell nextflow to run one alignment process at a time with maxForks 1.

    // NOTE 2: I need to decide on how to work with the possibility that users
    // will run multiple aligners (bwa-mem2, bwa, bsbolt, etc). We should still
    // only run one at a time.
    maxForks 1

    input:
    tuple val(sample_id), path(index_dir), val(index_prefix), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    shell:
    align = "bwa mem -t 10 -SP5M ${index_dir}/${index_prefix} ${fastq1} ${fastq2}"
    tobam = "samtools view -b -o ${sample_id}.bam"
    "${align} | ${tobam}"

    stub:
    "touch ${sample_id}.bam"
}

workflow Align {
    take:
    samples

    main:
    // Give options for BWA-MEM (for lower indexing requirements)
    // and BSBolt (for methylation + hi-c alignment)

    // Give option to modify number of threads to use for alignment
    // Give option for single-end reads

    samples
        | filter{it.datatype == "fastq"
                    && it.get("fastq1")
                    && it.get("fastq2")
        }
        | map{
            it.fastq1 = file(it.fastq1)
            it.fastq2 = file(it.fastq2)
            it
        }
        | set {fastq}
    
    fastq
        | branch {
            bwamem2: it.aligner == "bwa-mem2"
            bwamem: it.aligner == "bwa-mem"
            bsbolt: it.aligner == "bsbolt"
            error: true
        } | set{fastq}
    

    samples = transpack(
        BwaMem2Align,
        [fastq.bwamem2, samples],
        ["sample_id", "index_dir", "index_prefix", "fastq1", "fastq2"],
        ["sample_id", "sambam"],
        ["latest":"sambam"],
        "sample_id"
    )

    samples = transpack(
        BwaMemAlign,
        [fastq.bwamem, samples],
        ["sample_id", "index_dir", "index_prefix", "fastq1", "fastq2"],
        ["sample_id", "sambam"],
        ["latest":"sambam"],
        "sample_id"
    )


    if (params.general.get("last_step") == "align") {
        channel.empty() | set{samples}
    }

    emit:
    samples
}