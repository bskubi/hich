include {transpack; emptyOnLastStep} from './extraops.nf'

process ZcatHeadFastq {
    input:
    tuple val(id), path(fastq1), path(fastq2), val(n_reads)

    output:
    tuple val(id), path("${id}_R1.fastq.gz"), path("${id}_R2.fastq.gz")

    shell:
    lines = (n_reads as Integer) * 4

    "zcat ${fastq1} | head -n ${lines} | gzip -c > ${id}_R1.fastq.gz && zcat ${fastq2} | head -n ${lines} | gzip -c > ${id}_R2.fastq.gz"

    stub:
    "touch ${id}_R1.fastq.gz ${id}_R2.fastq.gz"
}

workflow HeadReads {
    take:
    samples

    main:
    // We should get a feature to downsample ingested bam and pairs files as well
    // to facilitate the --humid parameter
    
    samples
        | filter{it.datatype == "fastq" && it.get("n_reads") && it.get("n_reads").toString() != "all"}
        | set{head}

    samples = transpack(
        ZcatHeadFastq,
        [head, samples],
        ["id", "fastq1", "fastq2", "n_reads"],
        ["id", "fastq1", "fastq2"],
        [:],
        "id"
    )

    samples = emptyOnLastStep("headReads") ?: samples

    emit:
    samples
}