include {QCReads} from './qcHicReads.nf'
include {emptyOnLastStep; pack; skip} from './extraops.nf'

process PairtoolsFlipSort {
    publishDir params.general.publish.flip_sort ? params.general.publish.flip_sort : "results",
               saveAs: {params.general.publish.flip_sort ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::pairtools"
    container "bskubi/hich:latest"
    label 'doJobArray'
    label 'pairs'
    cpus 8

    input:
    tuple val(id), path(pairs), path(chromsizes), val(reshapeParams)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:

    reshapeParams = reshapeParams.join(" ")

    reshapeCmd = reshapeParams ? ["hich reshape ${reshapeParams}"] : []
    flipCmd = ["pairtools flip --chroms-path ${chromsizes}  --nproc-in ${task.cpus} --nproc-out ${task.cpus}"]
    sortCmd = ["pairtools sort --output ${id}.pairs.gz  --nproc-in ${task.cpus} --nproc-out ${task.cpus}"]

    if (reshapeParams) {
        reshapeCmd = reshapeCmd + ["--read_from ${pairs}"]
        reshapeCmd = [reshapeCmd.join(" ")]
    } else {
        flipCmd = flipCmd + ["${pairs}"]
        flipCmd = [flipCmd.join(" ")]
    }


    // Combine the individual commands, then join with a pipe to form the full command
    cmdParts = reshapeCmd + flipCmd + sortCmd
    cmd = cmdParts.join(" | ")

    // Execute the full command
    cmd

    stub:
    "touch ${id}.pairs.gz"
}

workflow IngestPairs {
    take:
        samples

    main:


    samples
        | filter{!skip("ingestPairs") && it.datatype == "pairs"}
        | map{tuple(it.id, it.pairs, it.chromsizes, it.reshapeParams)}
        | PairtoolsFlipSort
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{result}
    pack(samples, result) | set{samples}

    if ("ingestPairs" in params.general.get("qcAfter")) {
        samples = QCReads(samples, "ingestPairs")
    }

    samples = emptyOnLastStep("ingestPairs", samples)

    emit:
        samples
}