include {transpack} from './extraops.nf'

process BwaAlign {
    publishDir params.general.publish.bam ? params.general.publish.bam : "results",
               saveAs: {params.general.publish.bam ? it : null},
               mode: params.general.publish.mode

    conda "bwa bwa-mem2 samtools"
    container "bskubi/hich:latest"
    
    maxRetries 6
    memory {15.GB + 15.GB * (task.attempt-1)}

    // NOTE: Alignment speed is trivially parallelizeable and does not benefit
    // from running alignment in parallel multiple files at once. Each instance
    // of bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
    // tell nextflow to run one alignment process at a time with maxForks 1.

    // NOTE 2: I need to decide on how to work with the possibility that users
    // will run multiple aligners (bwa-mem2, bwa, bsbolt, etc). We should still
    // only run one at a time.
    maxForks 1

    input:
    tuple val(id), path(index_dir), val(index_prefix), path(fastq1), path(fastq2), val(aligner), val(threads), val(bwa_flags)

    output:
    tuple val(id), path("${id}.bam")

    shell:
    align = ""
    if (aligner in ["bwa-mem2", "bwa"]) {
        align = "${aligner} mem -t ${threads} ${bwa_flags} ${index_dir}/${index_prefix} ${fastq1} ${fastq2}"
    }
    
    tobam = "samtools view -b -o ${id}.bam"
    "${align} | ${tobam}"

    stub:
    "touch ${id}.bam"
}

workflow Align {
    take:
    samples

    main:
    // and BSBolt (for methylation + hi-c alignment)

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

    samples = transpack(
        BwaAlign,
        [fastq, samples],
        ["id", "index_dir", "index_prefix", "fastq1", "fastq2", "aligner", "aligner_threads", "bwa_flags"],
        ["id", "sambam"],
        ["latest":"sambam"],
        "id"
    )

    if (params.general.get("last_step") == "align") {
        channel.empty() | set{samples}
    }

    emit:
    samples
}