include {QCReads} from './qcHicReads.nf'
include {transpack; emptyOnLastStep; updateChannel} from './extraops.nf'

process PairtoolsParse2 {
    publishDir params.general.publish.parse ? params.general.publish.parse : "results",
               saveAs: {params.general.publish.parse ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::pairtools bioconda::samtools"
    container "bskubi/hich:latest"
    label 'doJobArray'
    cpus 8
    debug true

    input:
    tuple val(id), path(sambam), path(chromsizes), val(assembly), val(parseParams), val(reshapeParams)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:
    // The parseParams are typically a list of individual pairtools parse flags.
    // Join them separated by spaces to use in the parse2Cmd
    parseParams = parseParams.join(" ")
    reshapeParams = reshapeParams.join(" ")

    // Set up the individual commands in lists to make them easier to combine with pipes into a complete command
    // sambamba is both slower than samtools as of 2017, and also can't pipe to stdout, so we use samtools
    samSortCmd = ["samtools sort -n ${sambam}"]
    parse2Cmd = ["pairtools parse2 --assembly ${assembly} --chroms-path ${chromsizes} ${parseParams}"]
    reshapeCmd = reshapeParams ? ["hich reshape ${reshapeParams}"] : []
    pairsSortCmd = ["pairtools sort --output ${id}.pairs.gz --nproc-in ${task.cpus} --nproc-out ${task.cpus}"]

    // Combine the individual commands, then join with a pipe to form the full command
    cmdParts = samSortCmd + parse2Cmd + reshapeCmd + pairsSortCmd
    cmd = cmdParts.join(" | ")
    

    // Execute the full command
    cmd

    stub:
    "touch ${id}.pairs.gz"
}

workflow Parse {
    take:
    samples

    main:
    // This should give an error if the file does not exist
    
    // samples
    //     | filter{it.datatype == "sambam"}
    //     | map{
    //         sample ->
    //         if (!sample.sambam || !file(sample.sambam).exists()) {
    //             error "In sample with id ${sample.id}, sambam file is specified but does not exist"
    //         }

    //         sample
    //     }
    //     | set {sambam}

    samples | branch{yes: it.datatype in ["fastq", "sambam"]; no: true} | set {run}
    run.yes
        | map{tuple(it.id, it.sambam, it.chromsizes, it.assembly, it.parseParams, it.reshapeParams)}
        | PairtoolsParse2
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{runResult}
    updateChannel(run.yes, runResult) | concat(run.no) | set{samples}
    

    // samples = transpack(
    //     PairtoolsParse2,
    //     [sambam, samples],
    //     ["id", "sambam", "chromsizes", "assembly", "parseParams", "reshapeParams"],
    //     ["id", "pairs"],
    //     ["latest":"pairs"],
    //     "id",
    //     ["nullOk":["reshapeParams", "parseParams"]]
    //     )
    

    // It might be good to simplify these workflow control steps since they
    // are repeated frequently.
    if ("Parse" in params.general.get("qcAfter")) {
        QCReads(samples, "Parse")
    }

    samples = emptyOnLastStep("Parse", samples)

    emit:
        samples
}
