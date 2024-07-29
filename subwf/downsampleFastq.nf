include {transpack} from './extraops.nf'

process ZcatHeadFastq {
    input:
    tuple val(sample_id), path(fastq1), path(fastq2), val(n_reads)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz")

    shell:
    lines = n_reads * 4
    "zcat ${fastq1} | head -n ${lines} > ${sample_id}_R1.fastq.gz && zcat ${fastq2} | head -n ${lines} > ${sample_id}_R2.fastq.gz"
}

workflow HeadReads {
    take:
    samples

    main:
    samples
        | filter{it.datatype == "fastq" && it.get("n_reads") && it.get("n_reads").toString() != "all"}
        | map{print(it.get("fastq1").getClass()); it}
        | set{head}

    samples = transpack(
        ZcatHeadFastq,
        [head, samples],
        ["sample_id", "fastq1", "fastq2", "n_reads"],
        ["sample_id", "fastq1", "fastq2"],
        [:],
        "sample_id"
    )

    emit:
    samples
}