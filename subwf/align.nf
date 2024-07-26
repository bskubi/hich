include {JoinProcessResults} from './joinProcessResults.nf'
include {transpack; hashmapdiff} from './extraops.nf'

process BwaMem2Align {
    publishDir params.general.publish.bam ? params.general.publish.bam : "results",
               saveAs: {params.general.publish.bam ? it : null}

    container "bskubi/bwa-mem2"
    //maxRetries 4
    //memory {20.GB + 20.GB * task.attempt}

    // NOTE: Alignment speed is trivially parallelizeable and does not benefit
    // from running alignment in parallel multiple files at once. Each instance
    // of bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
    // tell nextflow to run one alignment process at a time with maxForks 1.
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

workflow jpr {
    take:
    fastq
    samples

    main:
    samples = JoinProcessResults(
        BwaMem2Align,
        [fastq, samples],
        ["sample_id", "index_dir", "index_prefix", "fastq1", "fastq2"],
        ["sample_id", "sambam"],
        ["sample_id"],
        false,
        "sambam")
    
    emit:
    samples
}

workflow tp {
    take:
    fastq
    samples

    main:
    samples = transpack(
        BwaMem2Align,
        [fastq, samples],
        ["sample_id", "index_dir", "index_prefix", "fastq1", "fastq2"],
        ["sample_id", "sambam"],
        ["latest":"sambam"],
        "sample_id"
    )

    emit:
    samples
}

workflow Align {
    take:
    samples

    main:
       
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
        BwaMem2Align,
        [fastq, samples],
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