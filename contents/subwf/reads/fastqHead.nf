include {emptyOnLastStep} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {withLog; stubLog} from '../util/logs.nf'

process ZcatHeadFastq {
    input:
    tuple val(id), path(fastq1), path(fastq2), val(n_reads)

    output:
    tuple val(id), path("${id}_R1.fastq.gz"), path("${id}_R2.fastq.gz")

    shell:
    lines = (n_reads as Integer) * 4

    cmd = "zcat '${fastq1}' | head -n ${lines} | gzip -c > '${id}_R1.fastq.gz' && zcat '${fastq2}' | head -n ${lines} | gzip -c > '${id}_R2.fastq.gz'"
    logMap = [task: "ZcatHeadFastq", input: [id: id, fastq1: fastq1, fastq2: fastq2], output: [fastq1: "${id}_R1.fastq.gz", fastq2: "${id}_R2.fastq.gz"]]
    withLog(cmd, logMap)

    stub:
    stub = "touch '${id}_R1.fastq.gz' '${id}_R2.fastq.gz'"
    lines = (n_reads as Integer) * 4

    cmd = "zcat '${fastq1}' | head -n ${lines} | gzip -c > '${id}_R1.fastq.gz' && zcat '${fastq2}' | head -n ${lines} | gzip -c > '${id}_R2.fastq.gz'"
    logMap = [task: "ZcatHeadFastq", input: [id: id, fastq1: fastq1, fastq2: fastq2], output: [fastq1: "${id}_R1.fastq.gz", fastq2: "${id}_R2.fastq.gz"]]
    stubLog(stub, cmd, logMap)
}

workflow FastqHead {
    take:
    samples

    main:
    // !TODO: Single end, plaintext, other compression, other seq data formats
    
    samples
        | filter{it.datatype == "fastq" && it.get("n_reads") && it.get("n_reads").toString() != "all" && it.get("fastq1") && it.get("fastq2")}
        | map{tuple(it.id, it.fastq1, it.fastq2, it.n_reads)}
        | ZcatHeadFastq
        | map{id, fastq1, fastq2 -> [id: id, fastq1: fastq1, fastq2: fastq2]}
        | set{result}

    keyUpdate(samples, result, "id") | set{samples}

    samples = emptyOnLastStep("fastqHead", samples)

    emit:
    samples
}