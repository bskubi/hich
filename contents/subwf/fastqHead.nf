include {emptyOnLastStep; pack} from './extraops.nf'

process ZcatHeadFastq {
    input:
    tuple val(id), path(fastq1), path(fastq2), val(n_reads)

    output:
    tuple val(id), path("${id}_R1.fastq.gz"), path("${id}_R2.fastq.gz")

    shell:
    lines = (n_reads as Integer) * 4

    "zcat '${fastq1}' | head -n ${lines} | gzip -c > '${id}_R1.fastq.gz' && zcat '${fastq2}' | head -n ${lines} | gzip -c > '${id}_R2.fastq.gz'"

    stub:
    "touch '${id}_R1.fastq.gz' '${id}_R2.fastq.gz'"
}

workflow FastqHead {
    take:
    samples

    main:
    // !TODO: Plaintext, other compression, other seq data formats
    
    samples
        | filter{it.datatype == "fastq" && it.get("n_reads") && it.get("n_reads").toString() != "all"}
        | map{tuple(it.id, it.fastq1, it.fastq2, it.n_reads)}
        | ZcatHeadFastq
        | map{id, fastq1, fastq2 -> [id: id, fastq1: fastq1, fastq2: fastq2]}
        | set{result}

    pack(samples, result) | set{samples}

    samples = emptyOnLastStep("fastqHead", samples)

    emit:
    samples
}