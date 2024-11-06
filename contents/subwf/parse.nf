include {QCReads} from './qcHicReads.nf'
include {emptyOnLastStep; pack; skip} from './extraops.nf'

process PairtoolsParse2 {
    publishDir params.general.publish.parse ? params.general.publish.parse : "results",
               saveAs: {params.general.publish.parse ? it : null},
               mode: params.general.publish.mode

    container params.general.hichContainer
    label 'doJobArray'
    label 'pairs'

    input:
    tuple val(id), path(sambam), path(chromsizes), val(assembly), val(parseParams), val(reshapeParams), val(flags)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:
    parseParams = parseParams ?: []
    reshapeParams = reshapeParams ?: []

    // Use minMapq as default, but override with manually specified --min-mapq
    if (flags.minMapq instanceof Integer && !parseParams.any{it.contains("--min-mapq")}) {
        parseParams += ["--min-mapq ${flags.minMapq}"]
    } 


    // The parseParams are typically a list of individual pairtools parse flags.
    // Join them separated by spaces to use in the parse2Cmd
    parseParams = parseParams.join(" ")
    reshapeParams = reshapeParams.join(" ")

    // Set up the individual commands in lists to make them easier to combine with pipes into a complete command
    // sambamba is both slower than samtools as of 2017, and also can't pipe to stdout, so we use samtools
    
    sortCmd = ["samtools sort -n '${sambam}'"]
    viewCmd = ["samtools view -b '${sambam}'"]

    parse2Cmd = ["pairtools parse2 --assembly '${assembly}' --chroms-path '${chromsizes}' ${parseParams}"]
    reshapeCmd = reshapeParams ? ["hich reshape ${reshapeParams}"] : []
    pairsSortCmd = ["pairtools sort --output '${id}.pairs.gz' --nproc-in ${task.cpus} --nproc-out ${task.cpus}"]

    // Combine the individual commands, then join with a pipe to form the full command
    parseCmd = parse2Cmd + reshapeCmd + pairsSortCmd
    sortParse = (sortCmd + parseCmd).join(" | ")
    viewParse = (viewCmd + parseCmd).join(" | ")

    
    cmd = "samtools view '${sambam}' | awk -F'\\t' '{ print \$1 }' | partitioned && ${viewParse} || ${sortParse}"

    // Execute the full command
    cmd

    stub:
    "touch '${id}.pairs.gz'"
}

workflow Parse {
    take:
    samples

    main:
    samples
        | filter{!skip("parse") && it.datatype in ["fastq", "sambam"]}
        | map{tuple(it.id, it.sambam, it.chromsizes, it.assembly, it.pairtoolsParse2Params, it.reshapeParams, it.subMap("minMapq"))}
        | PairtoolsParse2
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{result}
    pack(samples, result) | set{samples}

    // It might be good to simplify these workflow control steps since they
    // are repeated frequently.
    if ("parse" in params.general.get("qcAfter")) {
        QCReads(samples, "parse")
    }

    samples = emptyOnLastStep("parse", samples)

    emit:
        samples
}
